# https://github.com/casey/just

alias b := build
build:
  make solution

alias c := clean
clean:
  rm ./solution

[private]
test-core ID:
  cat ./test/open_tests/open_{{ID}}.txt > ./input.txt 
  ./solution < ./input.txt > ./output.txt
  diff ./output.txt ./test/open_tests/open_{{ID}}.output.txt -y | nl -w2 -s': ' | grep '|' && exit 1 || exit 0

alias t := test
test ID: build
  just test-core {{ID}}

alias tu := test-until
test-until ID: build
  #!/usr/bin/env bash
  for i in $(seq 1 {{ID}});
  do
    just test-core $i
  done
