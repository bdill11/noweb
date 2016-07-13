#!/bin/bash
cat filelist | tr -d '\\' | tr -d \\r > templist 
for f in `cat templist`
do
  if [[ "$f" =~ .*-.*\.m ]] 
  then mv "$f" `echo "$f" | sed -e s/-/_/g`
  fi
done 
rm templist
