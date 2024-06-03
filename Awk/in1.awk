BEGIN {
  if ( unit == "" ) {
    unit=5;
  }
}

/K-VECTORS FROM UNIT:*/ {
  l     = length($0)
  part1 = substr($0, 1, 20);
  part2 = substr($0, 22);
  lrest = l-22;
  printf "%20s%1d%"lrest"s\n", part1, unit, part2;
  exit;
}
/a*/ { print; }
