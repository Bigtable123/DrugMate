# Built-in model providers

`base.py` defines `ProviderProfile`: credentials, endpoints, model catalogs and
provider-specific request behavior. `builtin/` contains the profiles shipped
with this source release. `providers.__init__` loads a fixed list of these
modules lazily; external provider directories and package entry points are not
loaded.

The agent owns client creation, credential refresh, streaming, prompt caching
and conversation ordering. Profiles retain API-specific request fields and
message adaptations. Model selection and login use the existing
`get_provider_profile()` and `list_providers()` functions.

Use the built-in `custom` profile for an OpenAI-compatible endpoint. Keep keys
in the local application settings or process environment, outside the source
checkout.
