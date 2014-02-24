# ROI conversion -- old to new
#
# USAGE: awk -f roi_old_to_new <in.roi >out.roi

/^#/ { print $0 }
! /^#/ {

  printf ("%s\t%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $6, $5, $4, $7);
  i = 8
  for (; i <= NF - $5 + 1; i++)
  	printf("\t%s", $i);
  first = 0;
  for (; i <= NF; i++)
  {
    if (!first)
      printf(",%s", $i);
    else
      printf("%s", $i);
    first = 0;
  }
  printf("\n");
}
