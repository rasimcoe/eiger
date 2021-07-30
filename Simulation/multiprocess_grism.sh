for visit in 1 #can add 2 3 4 here
do 
	for module in a b
	do
		for quartile in 0 1 2 3
		do
			CMD="python simulate_grism_core.py $visit $module $quartile &"
			eval $CMD
		done
	done
done
 
