# https://github.com/casey/just

build:
  make solution

clean:
  rm ./solution

test ID: build
  ./solution < ./test/open_tests/open_{{ID}}.txt > ./output.txt
  diff ./output.txt ./test/open_tests/open_{{ID}}.output.txt -y | nl -w2 -s': ' | grep '|' && exit 1 || exit 0
