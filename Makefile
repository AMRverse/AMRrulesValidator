# Makefile for AMRrulesValidator

.PHONY: dev build clean test setup-resources install-with-resources

# Copy rules and install in editable mode
dev:
	@echo "📦 Installing package in editable mode..."
	pip install -e .
	@echo "📥 Downloading and setting up CARD resource files..."
	python -c "from amrrulesvalidator.utils.resources import ResourceManager; ResourceManager().download_card_archives()"

# Run setup-resources after installation
#install-with-resources: dev setup-resources

# Download and extract CARD resource files
#setup-resources:
#	@echo "📥 Downloading and setting up CARD resource files..."
#	python -c "from amrrulesvalidator.utils.resources import ResourceManager; ResourceManager().download_card_archives()"

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