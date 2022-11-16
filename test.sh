set -eu

g++ -O2 -o A A.cpp
cd tools
for i in $(seq 0 99)
do
  f=$(printf "%04d.txt" $i)
  cargo run --release --bin tester ../A < in/${f} > out/${f}
done
