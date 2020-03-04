export PATH="$PATH:~/bin"
echo $PATH

models=$HOME/DATA/listExclude
structures=$HOME/abymod/DATA/abpdblib
for mod in ${models}/*.pdb
do
    struc=${structures}/`basename $mod`
    rmsca=`profit -f cdrH3ca.pro $mod $struc | grep RMS`
    rmsall=`profit -f cdrH3all.pro $mod $struc | grep RMS`
    pdb=`basename $mod .pdb`
    echo "$pdb $rmsca $rmsall"
done
