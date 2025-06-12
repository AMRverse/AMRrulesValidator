# Makefile for AMRrulesValidator

.PHONY: dev build clean

# Copy rules and install in editable mode
dev:
	@echo "📦 Installing package in editable mode..."
	pip install -e .

# Build distribution packages (wheel + sdist)
build:
	@echo "🚀 Building package..."
	python -m build

# Clean generated artifacts
clean:
	@echo "🧹 Cleaning build artifacts..."
	rm -rf build dist *.egg-info