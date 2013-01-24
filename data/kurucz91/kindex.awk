# This awkscript splits the kurucz filenames into bits which give temperature
# gravity, etc.  It is used to create kurucz_index

awk '
{ 	i=index($1,"t");
	j=index($1,"g");
	tem=substr($1,i+1,j-i-1);
	i=index($1,"g");
	j=index($1,"k");
	gr=substr($1,i+1,j-i-1);
	printf "%s	%s	%s\n", $1,tem,gr
	}
' $1 >/tmp/file1
awk '{
	printf("kurucz91/%-30s %d          %3.1f\n",$1,$2,$3/10)
}
' /tmp/file1  >/tmp/file2
sort -g --key=2 tmp/file2 >kurucz91.ls
echo "Check kurucz91.ls for bad gravities"
rm /tmp/file1 /tmp/file2

