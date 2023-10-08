make && clear && cat ../data/stereo_l0_r9.raw | ./project 0 1 | aplay -c 2 -f S16_LE -r 48000
# python3 ../model/cppArrayAsPSD.py
python3 ../model/wavio.py 48000