#!/bin/bash

ENV_NAME='e3'

conda activate $ENV_NAME

pip freeze > requirements.txt

echo "Requirements updated"