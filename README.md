-- Install CMAKE >= 3.20

### Cloning with Submodules

To clone this repository and its submodules:

```
git clone --recurse-submodules https://github.com/AstridRivPar/Alignment-FMIndex.git
```

Or if you already cloned it:

```
git submodule update --init --recursive
```

## Building the project
**Create a build directory:**
  ```    
  mkdir build
  cd build
  cmake -DCMAKE_BUILD_TYPE=RELEASE ..
  ```
    
**To run the example**
  ```
  make run_FMIndex
  ```
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
