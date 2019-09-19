#!/bin/bash

grep '^##' M.honghuensis_1_blob.blobDB.table.txt;\
grep -v '^##' M.honghuensis_1_blob.blobDB.table.txt | \
column -t -s $'\t'


