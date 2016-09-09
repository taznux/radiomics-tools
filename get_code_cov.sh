
#!/bin/bash
for filename in `find . | egrep '\.cxx|\.cpp|\.hpp|\.hxx'`; 
do 
  gcov-5 -n -o . $filename > /dev/null; 
done
