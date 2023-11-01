#!/usr/bin/python3

from collections import Counter
import json

with open("Result/duplicates.txt", 'r') as F:
    All = F.read().replace('\t', ' ').split('\n')

Result = {}
[Result.update(dict(Counter([str(ii) + ": " +i.split('-')[0] for i in All[ii].split(' ')[1:]]))) for ii in range(len(All)-1)]

with open("Result/Dupli_count.json", "w") as f:
    json.dump(Result, f)
