
#!/bin/bash
for filename in `find . | egrep '\.cxx|\.cpp|\.hpp|\.hxx'`; 
do 
  gcov -n -o . $filename > /dev/null; 
done
