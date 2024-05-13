REGISTRY="ghcr.io/sdsc-ordes"
IMAGE_NAME="modo-api"
VERSION :=$(shell grep -E '^version += +' pyproject.toml | sed -E 's/.*= +//')

.PHONY: install
install: ## Install with the poetry and add pre-commit hooks
	@echo "🚀 Installing packages with poetry"
	@poetry install
	@poetry run pre-commit install

.PHONY: check
check: ## Run code quality tools.
	@echo "🚀 Checking Poetry lock file consistency with 'pyproject.toml': Running poetry lock --check"
	@poetry lock --check
	@echo "🚀 Linting code: Running pre-commit"
	@poetry run pre-commit run -a

.PHONY: doc
doc: ## Build sphinx documentation website locally
	@echo "📖 Building documentation"
	@cd docs
	@poetry run sphinx-build docs/ docs/_build

.PHONY: docker-build
docker-build: ## Build the modo-api client Docker image
	@echo "🐋 Building docker image"
	@docker build \
		--build-arg="VERSION_BUILD=$(VERSION)" \
		-t $(REGISTRY)/$(IMAGE_NAME):$(VERSION) .

.PHONY: test
test: ## Test the code with pytest
	@echo "🚀 Testing code: Running pytest"
	@poetry run pytest

.PHONY: help
help:
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

.DEFAULT_GOAL := help
