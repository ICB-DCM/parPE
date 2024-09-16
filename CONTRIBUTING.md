# How to contribute

We are happy about contributions to parPE in any form, be it new functionality,
documentation, bug reports, or anything else.

When planning to contribute to parPE, it may be best to first create an issue
at [https://github.com/ICB-DCM/parPE/issues](https://github.com/ICB-DCM/parPE/issues)
to outline plans to avoid any redundant work.

For any contribution, please create a pull request. By creating a pull request,
you agree on contributing your code/text/image/etc. under the license terms
stated in [https://github.com/ICB-DCM/parPE/blob/master/LICENSE](https://github.com/ICB-DCM/parPE/blob/master/LICENSE). 


## Style guide

For any code contributions, please follow the guidelines below.


### General

* All files and functions should come with file-level and function-level
  documentation.
  
* All new functionality should be covered by unit or integration tests. Runtime
  of those tests should be kept as short as possible. 


### C++

* We use C++14

* New contributions should follow the
  [Google C++ Style Guide](https://google.github.io/styleguide/cppguide.html).
  Existing code does not yet, but will be adapted accordingly.


### Python 

* We want to be compatible with Python 3.9

* For the Python code we want to follow 
  [PEP8](https://www.python.org/dev/peps/pep-0008/). Although this is not the
  case for all existing code, any new contributions should do so. 

* We use Python [type hints](https://docs.python.org/3/library/typing.html)
  for all functions.
