#!/usr/bin/env bash

bioawk -c fastx '{ print $name, length($seq) }' < $1

