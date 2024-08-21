# FoTruss

## Dataset

Download datasets and put them into the /Dataset folder.

Each multilayer graph should be stored in a .txt file, and the filename is considered the dataset name.

The first line of the .txt file should contain ``#layers, #vertex, #edges``.

Following the header, each row contains three integers, which are ``layer-id, node-id, node-id``.

An example file "homo" is provided, and new datasets should adhere to the same structure.

We also include a script to convert ``WikiTalk`` and ``StackOverflow`` from original snap format ([Wiki](https://snap.stanford.edu/data/wiki-Talk.html), [SO](https://snap.stanford.edu/data/sx-stackoverflow.html)) to multilayer graph. Refer to "ConvertWiki.py" and "CovertSO.py" for usage.

The two case study datasets are also included as ``dblp-ijcai-kdd-mod.txt`` and ``dblp-ijcai-kdd-mod-aaai.txt``. ``metainfo_id_author.txt`` shows a map from ``node-id`` to ``author name`` for ``dblp-ijcai-kdd-mod.txt``. For more information about how to convert DBLP dataset to multilayer graphs, check our previous code at [here](https://github.com/MDCGraph/DBLP-MLG).

The missing datasets are available at [FirmCore](https://github.com/joint-em/FTCS/tree/main/Code/Datasets) and [multilayer kCore](https://github.com/egalimberti/multilayer_core_decomposition).

## Usage

### Compile

FoTruss uses `gflags` `glog` and `fmt`. A compiler supporting C++17 features is required.

To build:

```shell
    mkdir build
    cd build
    cmake ..
    make
```

### Run

```shell
    ./build/FoTruss --routine [run/test] --dataset [dataset_name] --algo [algorithm_name]
```

## Example

### Decomposition

```shell
    ./build/FoTruss --routine run --dataset homo --algo decompsition
```

### Build Index

```shell
    ./build/FoTruss --routine run --dataset homo --algo index
```

### Query

```shell
    ./build/FoTruss --routine test --dataset RM --algo sota-gt
    ./build/FoTruss --routine test --dataset homo --algo sota-random
```