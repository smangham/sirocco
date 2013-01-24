# This awkscript splits the kurucz filenames into bits which give temperature
# gravity, etc.  It is used to create kurucz_index


{ 	i=index($1,"t");
	j=index($1,"g");
	tem=substr($1,i+1,j-i-1);
	i=index($1,"g");
	j=index($1,"k");
	gr=substr($1,i+1,j-i-1);
	printf "%s	%s	%s\n", $1,tem,gr
	}
