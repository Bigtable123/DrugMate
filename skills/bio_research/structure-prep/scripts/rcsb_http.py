"""Resilient HTTP helpers for RCSB requests used by structure preparation."""

from __future__ import annotations

import json
import logging
import socket
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from typing import Any


logger = logging.getLogger(__name__)

RCSB_REQUEST_TIMEOUT_SECONDS = 60
RCSB_MAX_RETRIES = 3
RCSB_MAX_ATTEMPTS = RCSB_MAX_RETRIES + 1
RCSB_RETRY_BACKOFF_SECONDS = (1.0, 2.0, 4.0)
RCSB_RETRYABLE_HTTP_CODES = {408, 429}


class RCSBRequestError(RuntimeError):
    """Describe an RCSB request failure after retry handling."""

    def __init__(
        self,
        *,
        url: str,
        attempts: int,
        max_retries: int,
        timeout_seconds: float,
        last_error: BaseException,
        status_code: int | None = None,
        retryable: bool,
    ) -> None:
        self.url = url
        self.attempts = attempts
        self.max_retries = max_retries
        self.timeout_seconds = timeout_seconds
        self.last_error_type = type(last_error).__name__
        self.last_error = str(last_error) or self.last_error_type
        self.status_code = status_code
        self.retryable = retryable
        self.retry_exhausted = retryable and attempts >= max_retries + 1
        super().__init__(
            f"RCSB request failed after {attempts} attempt(s) "
            f"(timeout={timeout_seconds:g}s, url={url}): "
            f"{self.last_error_type}: {self.last_error}"
        )

    def to_dict(self) -> dict[str, Any]:
        """Return JSON-serializable diagnostics for a skill failure payload."""
        return {
            "url": self.url,
            "attempts": self.attempts,
            "retries": max(0, self.attempts - 1),
            "max_retries": self.max_retries,
            "timeout_seconds": self.timeout_seconds,
            "status_code": self.status_code,
            "last_error_type": self.last_error_type,
            "last_error": self.last_error,
            "retryable": self.retryable,
            "retry_exhausted": self.retry_exhausted,
        }


class RCSBNotFoundError(RCSBRequestError):
    """Signal an HTTP 404 without retrying it."""


def request_bytes(
    request: urllib.request.Request | str,
    *,
    timeout_seconds: float = RCSB_REQUEST_TIMEOUT_SECONDS,
    max_retries: int = RCSB_MAX_RETRIES,
    sleep_fn: Callable[[float], None] | None = None,
) -> bytes:
    """Read one RCSB response with a 60-second timeout and up to 3 retries."""
    payload, _attempts = _request_bytes_with_attempts(
        request,
        timeout_seconds=timeout_seconds,
        max_retries=max_retries,
        sleep_fn=sleep_fn,
    )
    return payload


def request_json(
    request: urllib.request.Request | str,
    *,
    timeout_seconds: float = RCSB_REQUEST_TIMEOUT_SECONDS,
    max_retries: int = RCSB_MAX_RETRIES,
    sleep_fn: Callable[[float], None] | None = None,
) -> Any:
    """Read and decode one JSON response using the shared RCSB policy."""
    payload, attempts = _request_bytes_with_attempts(
        request,
        timeout_seconds=timeout_seconds,
        max_retries=max_retries,
        sleep_fn=sleep_fn,
    )
    url = _request_url(request)
    try:
        return json.loads(payload.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RCSBRequestError(
            url=url,
            attempts=attempts,
            max_retries=max_retries,
            timeout_seconds=timeout_seconds,
            last_error=exc,
            retryable=False,
        ) from exc


def _request_bytes_with_attempts(
    request: urllib.request.Request | str,
    *,
    timeout_seconds: float,
    max_retries: int,
    sleep_fn: Callable[[float], None] | None,
) -> tuple[bytes, int]:
    if timeout_seconds <= 0:
        raise ValueError("timeout_seconds must be greater than zero")
    if max_retries < 0:
        raise ValueError("max_retries must be non-negative")

    url = _request_url(request)
    sleeper = sleep_fn or time.sleep
    max_attempts = max_retries + 1

    for attempt in range(1, max_attempts + 1):
        try:
            with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
                return response.read(), attempt
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                raise RCSBNotFoundError(
                    url=url,
                    attempts=attempt,
                    max_retries=max_retries,
                    timeout_seconds=timeout_seconds,
                    last_error=exc,
                    status_code=exc.code,
                    retryable=False,
                ) from exc
            retryable = (
                exc.code in RCSB_RETRYABLE_HTTP_CODES
                or 500 <= exc.code <= 599
            )
            if not retryable:
                raise RCSBRequestError(
                    url=url,
                    attempts=attempt,
                    max_retries=max_retries,
                    timeout_seconds=timeout_seconds,
                    last_error=exc,
                    status_code=exc.code,
                    retryable=False,
                ) from exc
            last_error: BaseException = exc
            status_code: int | None = exc.code
        except (TimeoutError, socket.timeout, urllib.error.URLError) as exc:
            last_error = exc
            status_code = None

        if attempt >= max_attempts:
            raise RCSBRequestError(
                url=url,
                attempts=attempt,
                max_retries=max_retries,
                timeout_seconds=timeout_seconds,
                last_error=last_error,
                status_code=status_code,
                retryable=True,
            ) from last_error

        delay = _retry_delay(attempt)
        logger.warning(
            "RCSB request attempt %d/%d failed; retrying in %.1fs: %s: %s",
            attempt,
            max_attempts,
            delay,
            type(last_error).__name__,
            last_error,
        )
        sleeper(delay)

    raise AssertionError("RCSB retry loop exited unexpectedly")


def _request_url(request: urllib.request.Request | str) -> str:
    if isinstance(request, urllib.request.Request):
        return request.full_url
    return str(request)


def _retry_delay(failed_attempt: int) -> float:
    index = min(failed_attempt - 1, len(RCSB_RETRY_BACKOFF_SECONDS) - 1)
    return RCSB_RETRY_BACKOFF_SECONDS[index]
