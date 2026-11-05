#!/usr/bin/env python3
"""Create bounded campaigns and bind Agent-authored decisions to runnable configs."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
from campaign import initialize_campaign, authorize_campaign, confirm_final, abandon_run, bind_round_decision, _atomic_json, _load


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    init = commands.add_parser("init")
    init.add_argument("--state", type=Path, required=True)
    init.add_argument("--autonomy", choices=("A0", "A1", "A2", "A3"), default="A1")
    init.add_argument("--authorized-through", choices=("R1", "R2", "R3"))
    init.add_argument("--max-runs", type=int, default=6)
    init.add_argument("--authorization-ref", help="Reference to actual prior user authorization; omission records a proposed intent only")
    authorize = commands.add_parser("authorize")
    authorize.add_argument("--state", type=Path, required=True)
    authorize.add_argument("--authorization-ref", required=True)
    authorize.add_argument("--autonomy", choices=("A0", "A1", "A2", "A3"))
    authorize.add_argument("--authorized-through", choices=("R1", "R2", "R3"))
    authorize.add_argument("--max-runs", type=int)
    confirm = commands.add_parser("confirm-final")
    confirm.add_argument("--state", type=Path, required=True)
    confirm.add_argument("--run-id", required=True)
    confirm.add_argument("--authorization-ref", required=True)
    abandon = commands.add_parser("abandon")
    abandon.add_argument("--state", type=Path, required=True)
    abandon.add_argument("--run-id", required=True)
    abandon.add_argument("--reason", required=True)
    plan = commands.add_parser("plan")
    plan.add_argument("--config", type=Path, required=True)
    plan.add_argument("--decision", type=Path, required=True)
    plan.add_argument("--output", type=Path, required=True)
    status = commands.add_parser("status")
    status.add_argument("--state", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "init":
            result = initialize_campaign(args.state, autonomy=args.autonomy,
                authorized_through=args.authorized_through, max_runs=args.max_runs,
                authorization_ref=args.authorization_ref)
        elif args.command == "authorize":
            result = authorize_campaign(args.state, autonomy=args.autonomy,
                authorized_through=args.authorized_through, max_runs=args.max_runs,
                authorization_ref=args.authorization_ref)
        elif args.command == "confirm-final":
            result = confirm_final(args.state, run_id=args.run_id,
                authorization_ref=args.authorization_ref)
        elif args.command == "abandon":
            result = abandon_run(args.state, run_id=args.run_id, reason=args.reason)
        elif args.command == "plan":
            config = _load(args.config)
            output = args.output.expanduser().resolve()
            expected = (config.get("project") or {}).get("round_decision_file")
            if not expected or output != Path(expected).expanduser().resolve():
                raise ValueError("--output must equal config.project.round_decision_file")
            result = bind_round_decision(config, _load(args.decision))
            _atomic_json(output, result)
        else:
            result = _load(args.state)
        print(json.dumps(result, indent=2, ensure_ascii=False))
        return 0
    except (OSError, ValueError, TypeError, KeyError) as exc:
        print(json.dumps({"success": False, "error": str(exc)}, ensure_ascii=False), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
