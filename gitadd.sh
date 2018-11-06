#!/bin/bash
branch="dav_kayamash"
file1="Newrun.C"
file2="NewEventDisplay.C"
file3="NewEventDisplay.h"
file4="HowToUse.txt"
file5="gitadd.sh"
add="git add "
message="add list"
push="git push origin "

eval $add$file1
eval $add$file2
eval $add$file3
eval $add$file4
eval $add$file5
git commit -m "${message}"
eval $push$branch
