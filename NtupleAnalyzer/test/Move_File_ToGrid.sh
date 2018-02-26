
if [ "$1" == "" ]
then
echo "Missing arg : file name"
elif [ "$2" == "" ]
then
echo "Missing arg : destination"
fi

#Copy file from scratch1 to grid
gfal-mkdir srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleAnalyzer/$2
gfal-copy -r $1 srm://sbgse1.in2p3.fr:8446/dpm/in2p3.fr/home/cms/phedex/store/user/ntonon/NtupleAnalyzer/$2/.

#Then remove file from scratch1
rm $1
