#Script to automatically copy the filelists for jobs which have failed, so that they can then be run another time
#Based on a grep (e.g. "crash") to match jobs which failed

#echo "--- Usage : ./Copy_FileLists_Failed.sh [logdirpath] [word to match for failed jobs]"

dir_allLists="./lists"
dir_listsToRun="./lists_tHq"

#Check datacard argument
if [[ $1 == "" || $2 == "" ]] ; then
    echo "--- Usage : ./Copy_FileLists_Failed.sh [logdirpath] [word to match for failed jobs]"
    exit
fi

#cut -d'X' -f1 : cuts in 2 parts at the 'X' delimiter, and returns first part
grep -irl "$2" $1 | while read line
do
  #echo ${line}
  sampStrip=$(echo ${line} | cut -d'/' -f2 | cut -d'.' -f1)
  #echo ${sampStrip}
  file="$dir_allLists/$sampStrip.txt"

  echo ""
  echo "cp $file $dir_listsToRun/."
  cp $file $dir_listsToRun/.
done  
