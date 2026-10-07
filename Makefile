.PHONY: install validate test audit

install:
	python -m pip install -e '.[enrichment,notebook,test]'

validate:
	sma-xspecies validate

test:
	pytest

audit: validate test
	sma-xspecies report

