from pathlib import Path
import tomllib

## Define repo file structure
working_dir = Path(__file__).parent.resolve()
# Data
data_dir = working_dir / "data"
raw_data_dir = data_dir / "raw"
processed_data_dir = data_dir / "processed"
# Figures
figure_dir = working_dir / "figures"

## Load user-defined vars from TOML file
try:
    with (working_dir / "config.toml").open("rb") as fp:
        user_config = tomllib.load(fp)
except FileNotFoundError:
    raise RuntimeError("Config file not found; copy `config.toml.example` to `config.toml` and fill in user-specific paths.")