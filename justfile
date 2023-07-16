# https://github.com/casey/just

alias b := build
build:
  make solution

alias c := clean
clean:
  rm -f ./solution
  rm -f massif.txt
  rm -f massif.out.*
  rm -f callgrind.out.*
  rm -f input.txt
  rm -f output.txt
  touch input.txt
  touch output.txt

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

valgrind ID: build
  valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --log-fd=9 9>&1 ./solution < ./test/open_tests/open_{{ID}}.txt 1>./output.txt

massif ID: build
  cat ./test/open_tests/open_{{ID}}.txt > ./input.txt
  valgrind --tool=massif --log-fd=9 9>&1 ./solution < ./input.txt 1>./output.txt
  ms_print massif.out.* > massif.txt
  rm massif.out.*

callgrind ID: build
  cat ./test/open_tests/open_{{ID}}.txt > ./input.txt
  valgrind --tool=callgrind --log-fd=9 9>&1 ./solution < ./input.txt 1>./output.txt
  kcachegrind callgrind.out.*
  rm callgrind.out.*
