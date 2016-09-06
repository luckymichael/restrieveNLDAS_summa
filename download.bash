
for h in {1..23}
do
	hh="0$h"
	echo ${hh: -2}
done

# slightly malformed input data
input_start='2016-02-20'
input_end='2016-03-02'

# After this, startdate and enddate will be valid ISO 8601 dates,
# or the script will have aborted when it encountered unparseable data
# such as input_end=abcd
startdate=$(date -I -d "$input_start") || exit -1
enddate=$(date -I -d "$input_end")     || exit -1

d="$startdate"
while [ "$d" != "$enddate" ]; do 
  echo $d
  jd=`date --date=$d +%j` # julian day
  yy=`date --date=$d +%Y` # year
  mm=`date --date=$d +%m` # month
  dd=`date --date=$d +%d` # day

  echo $jd $yy $mm $dd
  
  # increase one day
  d=$(date -I -d "$d + 1 day")
 
done
