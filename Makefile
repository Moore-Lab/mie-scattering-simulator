# mieinfo — build/validate targets (INTERFACES.md §10). `make validate` is the
# integration gate the orchestrator runs before every merge to develop.
.PHONY: install test validate optimize compare report lint

install:
	pip install -e ".[dev]"

test:               # fast PR lane
	pytest -q -m "not slow"

validate:           # full validation gates; nonzero exit on failure
	python -m mieinfo.cli validate

lint:
	ruff check mieinfo && mypy mieinfo

optimize:
	python -m mieinfo.cli optimize $(CONFIG)

compare:
	python -m mieinfo.cli compare

report:
	python -m mieinfo.cli report
