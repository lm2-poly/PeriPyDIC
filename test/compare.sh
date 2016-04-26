md5_stable=($(md5sum $1))
md5_new=($(md5sum $1))

if [ "$md5_stable" == "$md5_new" ]; then
		echo "Test passed"
else
		echo "test failed"
fi
