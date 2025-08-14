/// Takes 2 segmentations, finds the best 1-1 segment correspondence,
/// renumbers the segments in the second segmentation, and saves it

#include <stdio.h>
#include <tchar.h>

#include <iostream>
#include <map>
#include <set>


using namespace std;

#define max(a, b) (((a) < (b)) ? (b) : (a))
#define min(a, b) (((a) < (b)) ? (a) : (b))
#define minmax(a, b, c) (min(b, max(a, c)))
#define sqr(x) ((x) * (x))

double euclid(int a, int b) { return sqrt(sqr((double)a) + sqr((double)b)); }

typedef uint16_t pixtype;

typedef struct {
    int hor, ver;
} dmpixtype;  /// type for distance map with horizontal and vertical distance separated for fast evaluation of euclid metric

// connects overlapping segments from different segmentations
typedef std::set<std::pair<int, int>> bind_t;

// Clamps variable and saves the result back.
// Does not check if `variable`, `min`, or `max` are not NAN.
// \returns false if `min` > `max` or does not change `variable`
template <typename T1, typename T2, typename T3>
bool stdClampInPlace(T1& variable, const T2& min, const T3& max) noexcept {
    if (min > max) return false;
    if (variable < min) {
        variable = min;
        return true;
    }
    if (variable > max) {
        variable = max;
        return true;
    }
    return false;
}

void export_image(char* filename, pixtype* data, int width, int height) {
    FILE* file = _fsopen(filename, "wb", SH_DENYWR);
    fwrite(data, sizeof(pixtype), width * height, file);
    fclose(file);
}

void findGrouping121(pixtype* in1, pixtype* in2, int width, int height, bind_t& bind) {
    typedef std::pair<int, int> segm_pair_t;
    typedef std::map<segm_pair_t, int> intersect_t;

    /// find sizes of intersections
    intersect_t inter;
    {
        segm_pair_t pair;
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                pair.first = in1[j * width + i];
                pair.second = in2[j * width + i];
                ++inter[pair];
            }
        }
    }

    std::set<int> processed1, processed2;

    /// find bindings
    while (true) {
        /// find max intersection
        int max = 0;
        segm_pair_t max_pair;
        for (const auto& it : inter) {
            it.second;
            if (it.second <= max) continue;
            if (processed1.find(it.first.first) != processed1.end()) continue;
            if (processed2.find(it.first.second) != processed2.end()) continue;

            max_pair = it.first;
            max = it.second;
        }

        if (max == 0) break;  /// no more intersections
        /// it's a new maximal intersection
        bind.insert(max_pair);
        processed1.insert(max_pair.first);
        processed2.insert(max_pair.second);
    }
}

/// uses segment numbers from the 1st segmentation to "rewrite" segment numbers in the 2nd segmentation
/// \param in input segmentation, is not changed
void rewriteSegments(const pixtype* in, pixtype* out, int width, int height, const bind_t& bind) {
    std::map<int, int> map;  ///< maps original segment numbers to new ones
    int last_new = 0;        ///< last new segment number

    /// initialize map
    for (auto const& it : bind) map[it.second] = it.first;
    map[0] = 0;  ///< keep background

    for (int i = 0; i < width * height; i++) {
        if (map.find(in[i]) != map.end()) {
            out[i] = map[in[i]];
            continue;
        }

        /// find unused segment number
        for (int j = last_new + 1; j < INT_MAX; j++) {
            if (map.find(j) != map.end()) continue;
            last_new = j;
            break;
        }

        /// save and use new segment number
        map[in[i]] = last_new;
        out[i] = last_new;
    }
}

bool readImages(int width, int height, char* file1, char* file2, pixtype* segm1, pixtype* segm2) {
    FILE* fileseg;

    /// read 1st file
    fileseg = _fsopen(file1, "rb", SH_DENYWR);
    if (fileseg == 0) {
        cout << "Cannot open " << file1;
        return false;
    }

    fread(segm1, sizeof(pixtype), width * height, fileseg);
    fclose(fileseg);

    /// read 2nd file
    fileseg = _fsopen(file2, "rb", SH_DENYWR);
    if (fileseg == 0) {
        cout << "Cannot open " << file2;
        return false;
    }

    fread(segm2, sizeof(pixtype), width * height, fileseg);
    fclose(fileseg);
    return true;
}

void processSegmentations(int width, int height, char* file1, char* file2, char* file3) {
    pixtype* segm1 = (pixtype*)malloc(width * height * sizeof(pixtype));
    pixtype* segm2 = (pixtype*)malloc(width * height * sizeof(pixtype));
    pixtype* out = (pixtype*)malloc(width * height * sizeof(pixtype));

    if (!readImages(width, height, file1, file2, segm1, segm2)) return;

    bind_t bind;
    findGrouping121(segm1, segm2, width, height, bind);
    rewriteSegments(segm2, out, width, height, bind);
    export_image(file3, out, width, height);
}

int main(int argc, char* argv[]) {
    if (argc != 6) cout << "Provide width, height, 2 paths to files with segmentations, and one for output (RAW, 16 bit).";
    int width = atoi(argv[1]);
    int height = atoi(argv[2]);
    if (width <= 0) cout << "Invalid width.";
    if (height <= 0) cout << "Invalid height.";

    processSegmentations(width, height, argv[3], argv[4], argv[5]);
    return 0;
}