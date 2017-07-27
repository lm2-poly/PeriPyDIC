## Tests
cd ../test
echo "--1D direction x+"
python ../pd_dic.py -i input_elas_1D_x+.yaml -t pd > 1D+.dat
sed -i '$ d' 1D+.dat
DIFF=$(diff 1D.res 1D+.dat )
if [ "$DIFF" != "" ]
then
        echo "Test failed"
else
        echo "Test passed"
fi
echo "--1D direction x-"
python ../pd_dic.py -i input_elas_1D_x-.yaml -t pd > 1D-.dat
sed -i '$ d' 1D-.dat
DIFF=$(diff 1D.res 1D-.dat )
if [ "$DIFF" != "" ]
then
        echo "Test failed"
else
        echo "Test passed"
fi

