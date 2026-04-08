#!/bin/bash
set -o errexit

ruff check golfy/ tests/

echo 'Passes ruff check'
