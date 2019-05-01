#!/usr/bin/env bash

# Note: 
# This is just a template for genral use cases @TreX, it does not check for trueness of arguments. 
# User can supply at most 2 non acceptable values at atleast $1 $2 or $3 to execute true arguments. 

set -e

display_usage(){
  echo "------------------------------------------------------------------------------------------------------------------"
  echo "run the script using the following syntax:"
  echo "    bash" $0 "<-k1, -k2 > <-k3, -k4 > <-k5, -k6 >"
  echo ""
  echo "------------------------------------------------------------------------------------------------------------------"
  exit 1
}

if [[ $1 = "--help" ]] || [[ $1 = "-h" ]]; then
  #statements
  display_usage
  exit 1
fi

# check for missing arguments
if [[ -z $1 ]] || [[ -z $2 ]] || [[ -z $3 ]]; then
  #statements
  echo "------------------------------------------------------------------------------------------------------------------"
  echo "All \$1 \$2 & \$3 arguments needed"
  display_usage
  echo "------------------------------------------------------------------------------------------------------------------"
  exit 1
fi



func1(){
  echo "this is func1"
}
func2(){
  echo "this is func2"
}
func3(){
  echo "this is func3"
}
func4(){
  echo "this is func4"
}
func5(){
  echo "this is func5"
}
func6(){
  echo "this is func6"
}

raise_error() {
  echo "-------------------------------------------------------------------"
  local error_message="$@"
  echo "${error_message}" 1>&2;
  echo "-------------------------------------------------------------------"
}

case $1 in
    -k1|--knit1)
      func1
      ;;
    -k2|--knit2)
      func2
      ;;
    *)
      raise_error "Unknown/Missing argument(s): ${1}"
      exit 1
      ;;
esac

case $2 in
    -k3|--knit3)
      func3
      ;;
    -k4|--knit4)
      func4
      ;;
    *)
      raise_error "Unknown/Missing argument(s) at \$2 : ${2}"
      exit 1
      ;;
esac

case $3 in
    -k5|--knit5)
      func5
      ;;
    -k6|--knit6)
      func6
      ;;
    *)
      raise_error "Unknown/Missing argument(s) at \$3 : ${3}"
      exit 1
      ;;
esac
