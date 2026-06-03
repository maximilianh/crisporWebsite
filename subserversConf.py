"""Subserver registry — single source of truth shared by startSubServers.sh
and crispor.py. Each entry maps a module name (file <name>.py exposing run(data))
to its port and the venv directory whose python3 binary should host it."""

SUBSERVERS = {
    "runDeepBe":     {"port": 8001, "venv": "venvDeepBe"},
    "runForecastBe": {"port": 8002, "venv": "venvForecastBe"},
    # "runCrisprOnBe": {"port": 8003, "venv": "venvCrisprOnBe"}
}
