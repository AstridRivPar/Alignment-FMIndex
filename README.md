-- Requires CMAKE >= 3.20

### Cloning with Submodules

To clone this repository and its submodules:

```
git clone --recurse-submodules https://github.com/AstridRivPar/Alignment-FMIndex.git
```

If you've already cloned it without submodules, initialize submodules manually:

```
git submodule update --init --recursive
```

## Building the project
**Create a build directory:**
  ```    
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=RELEASE ..
  make
  ```
    
**To run the example**
Before running the example, create a directory named OutputFiles in the root of the repository (i.e., at the same level as the build/ folder):
```
mkdir OutputFiles
```
Then from the build directory, run:
  ```
  make run_FMIndex
  ```
The example will write its results to the OutputFiles folder.


**Parameters**
Optional:
   - `-p`: (Optional) Returns partial matches. **Default: return complete traces**.
  - `-f`: (Optional) Direction of the search. **Default: search backward**.
Required:
  - `-i <index_file>` - File with the runs in the model as strings. Each run must be in a separate line.
  - `-q <queries_file>` - File with the unique traces in the log as strings. Each trace is followed by the number of its instances in the log (Check example). 
  - `-o <output_file>` - File to output the alignments.

**Output**
The output of the code is in the following format:

`<trace>, <number of instances>, <number of moves>, <number of optimal alignments found>, <list of alignments>`
