# Makefile for AMRrulesValidator

.PHONY: dev build clean test setup-resources install-with-resources

# Copy rules and install in editable mode
dev:
	@echo "📦 Installing package in editable mode..."
	pip install -e .
	@echo "📥 Downloading and setting up resource files..."
	python -c "from amrrulesvalidator.utils.resources import ResourceManager; ResourceManager().setup_all_resources()"
	@echo "✓ Setup complete, resources include CARD ontology, drug information, and AMRFinderPlus data"

# Run tests
test:
	@echo "🧪 Running tests..."
	pytest -v

# Build distribution packages (wheel + sdist)
build:
	@echo "🚀 Building package..."
	python -m build

# Clean generated artifacts
clean:
	@echo "🧹 Cleaning build artifacts..."
	rm -rf build dist *.egg-info

# Setup resources explicitly
setup-resources:
	@echo "📥 Downloading and setting up resource files..."
	python -c "from amrrulesvalidator.utils.resources import ResourceManager; rm = ResourceManager(); success = rm.setup_all_resources(); exit(0 if success else 1)"

# Install with resources
install-with-resources: dev