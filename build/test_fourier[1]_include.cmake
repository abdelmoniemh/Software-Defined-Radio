if(EXISTS "/home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build/test_fourier[1]_tests.cmake")
  include("/home/moniem/code/3DY4/3dy4-project-group11-prj-wednesday/build/test_fourier[1]_tests.cmake")
else()
  add_test(test_fourier_NOT_BUILT test_fourier_NOT_BUILT)
endif()
