#include <stdio.h>
#include <tchar.h>

#include <iostream>

#include "stdafx.h"

// #define EXPORT_DISTANCE_MAP

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

typedef struct {
    int x1, y1, x2, y2;
} dmroi;

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

void export_distance_map(char* filename, dmpixtype* data, int width, int height) {
    pixtype* ex = (pixtype*)malloc(height * width * sizeof(pixtype));

    for (int j = 0; j < height; ++j) {
        for (int i = 0; i < width; ++i) {
            ex[j * width + i] = pixtype(3 * euclid(data[j * width + i].hor, data[j * width + i].ver));
        }
    }
    export_image(filename, ex, width, height);
}

/// computes distance map from pre-filled data dm; computation done only in region of interest
void computeDistanceMapROI(dmpixtype* dm, int height, int width, dmroi roi) {
    /// check for minimal and maximal values of ROI
    stdClampInPlace(roi.x1, 0, width - 1);
    stdClampInPlace(roi.x2, 0, width - 1);
    stdClampInPlace(roi.y1, 0, height - 1);
    stdClampInPlace(roi.y2, 0, height - 1);

    /// compute distance map - top down
    for (int j = roi.x1 + 1; j <= roi.x2; j++) {
        if ((sqr((uint64_t)dm[roi.y1 * width + j - 1].hor + 1) + sqr((uint64_t)dm[roi.y1 * width + j - 1].ver)) <
            (sqr((uint64_t)dm[roi.y1 * width + j].hor) + sqr((uint64_t)dm[roi.y1 * width + j].ver))) {
            dm[roi.y1 * width + j].hor = dm[roi.y1 * width + j - 1].hor + 1;
            dm[roi.y1 * width + j].ver = dm[roi.y1 * width + j - 1].ver;
        }
    }
    for (int i = roi.y1 + 1; i <= roi.y2; i++) {
        if ((sqr((uint64_t)dm[(i - 1) * width + roi.x1].hor) + sqr((uint64_t)dm[(i - 1) * width + roi.x1].ver + 1)) <
            (sqr((uint64_t)dm[i * width + roi.x1].hor) + sqr((uint64_t)dm[i * width + roi.x1].ver))) {
            dm[i * width + roi.x1].hor = dm[(i - 1) * width + roi.x1].hor;
            dm[i * width + roi.x1].ver = dm[(i - 1) * width + roi.x1].ver + 1;
        }
        for (int j = roi.x1 + 1; j <= roi.x2; j++) {
            if ((sqr((uint64_t)dm[(i - 1) * width + j].hor) + sqr((uint64_t)dm[(i - 1) * width + j].ver + 1)) <
                (sqr((uint64_t)dm[i * width + j].hor) + sqr((uint64_t)dm[i * width + j].ver))) {
                dm[i * width + j].hor = dm[(i - 1) * width + j].hor;
                dm[i * width + j].ver = dm[(i - 1) * width + j].ver + 1;
            }
            if ((sqr((uint64_t)dm[i * width + j - 1].hor + 1) + sqr((uint64_t)dm[i * width + j - 1].ver)) <
                (sqr((uint64_t)dm[i * width + j].hor) + sqr((uint64_t)dm[i * width + j].ver))) {
                dm[i * width + j].hor = dm[i * width + j - 1].hor + 1;
                dm[i * width + j].ver = dm[i * width + j - 1].ver;
            }
        }
    }

    /// compute distance map - bottom up
    for (int j = roi.x2 - 1; j >= roi.x1; j--) {
        if ((sqr((uint64_t)dm[roi.y2 * width + j + 1].hor + 1) + sqr((uint64_t)dm[roi.y2 * width + j + 1].ver)) <
            (sqr((uint64_t)dm[roi.y2 * width + j].hor) + sqr((uint64_t)dm[roi.y2 * width + j].ver))) {
            dm[roi.y2 * width + j].hor = dm[roi.y2 * width + j + 1].hor + 1;
            dm[roi.y2 * width + j].ver = dm[roi.y2 * width + j + 1].ver;
        }
    }
    for (int i = roi.y2 - 1; i >= roi.y1; i--) {
        if ((sqr((uint64_t)dm[(i + 1) * width + roi.x2].hor) + sqr((uint64_t)dm[(i + 1) * width + roi.x2].ver + 1)) <
            (sqr((uint64_t)dm[i * width + roi.x2].hor) + sqr((uint64_t)dm[i * width + roi.x2].ver))) {
            dm[i * width + roi.x2].hor = dm[(i + 1) * width + roi.x2].hor;
            dm[i * width + roi.x2].ver = dm[(i + 1) * width + roi.x2].ver + 1;
        }
        for (int j = roi.x2 - 1; j >= roi.x1; j--) {
            if ((sqr((uint64_t)dm[(i + 1) * width + j].hor) + sqr((uint64_t)dm[(i + 1) * width + j].ver + 1)) <
                (sqr((uint64_t)dm[i * width + j].hor) + sqr((uint64_t)dm[i * width + j].ver))) {
                dm[i * width + j].hor = dm[(i + 1) * width + j].hor;
                dm[i * width + j].ver = dm[(i + 1) * width + j].ver + 1;
            }
            if ((sqr((uint64_t)dm[i * width + j + 1].hor + 1) + sqr((uint64_t)dm[i * width + j + 1].ver)) <
                (sqr((uint64_t)dm[i * width + j].hor) + sqr((uint64_t)dm[i * width + j].ver))) {
                dm[i * width + j].hor = dm[i * width + j + 1].hor + 1;
                dm[i * width + j].ver = dm[i * width + j + 1].ver;
            }
        }
    }

    /// compute distance map - top down
    for (int j = roi.x1 + 1; j <= roi.x2; j++) {
        if ((sqr((uint64_t)dm[roi.y1 * width + j - 1].hor + 1) + sqr((uint64_t)dm[roi.y1 * width + j - 1].ver)) <
            (sqr((uint64_t)dm[roi.y1 * width + j].hor) + sqr((uint64_t)dm[roi.y1 * width + j].ver))) {
            dm[roi.y1 * width + j].hor = dm[roi.y1 * width + j - 1].hor + 1;
            dm[roi.y1 * width + j].ver = dm[roi.y1 * width + j - 1].ver;
        }
    }
    for (int i = roi.y1 + 1; i <= roi.y2; i++) {
        if ((sqr((uint64_t)dm[(i - 1) * width + roi.x1].hor) + sqr((uint64_t)dm[(i - 1) * width + roi.x1].ver + 1)) <
            (sqr((uint64_t)dm[i * width + roi.x1].hor) + sqr((uint64_t)dm[i * width + roi.x1].ver))) {
            dm[i * width + roi.x1].hor = dm[(i - 1) * width + roi.x1].hor;
            dm[i * width + roi.x1].ver = dm[(i - 1) * width + roi.x1].ver + 1;
        }
        for (int j = roi.x1 + 1; j <= roi.x2; j++) {
            if ((sqr((uint64_t)dm[(i - 1) * width + j].hor) + sqr((uint64_t)dm[(i - 1) * width + j].ver + 1)) <
                (sqr((uint64_t)dm[i * width + j].hor) + sqr((uint64_t)dm[i * width + j].ver))) {
                dm[i * width + j].hor = dm[(i - 1) * width + j].hor;
                dm[i * width + j].ver = dm[(i - 1) * width + j].ver + 1;
            }
            if ((sqr((uint64_t)dm[i * width + j - 1].hor + 1) + sqr((uint64_t)dm[i * width + j - 1].ver)) <
                (sqr((uint64_t)dm[i * width + j].hor) + sqr((uint64_t)dm[i * width + j].ver))) {
                dm[i * width + j].hor = dm[i * width + j - 1].hor + 1;
                dm[i * width + j].ver = dm[i * width + j - 1].ver;
            }
        }
    }
}

/// set distance of border pixels to zero, other to infinity
void createBorderDistanceMapROI(pixtype* im, dmpixtype* dm, int height, int width, pixtype seg, dmroi roi) {
    /// check for minimal and maximal values of ROI
    stdClampInPlace(roi.x1, 0, width - 1);
    stdClampInPlace(roi.x2, 0, width - 1);
    stdClampInPlace(roi.y1, 0, height - 1);
    stdClampInPlace(roi.y2, 0, height - 1);

    /// set maximal distances
    for (int i = roi.y1; i <= roi.y2; i++) {
        for (int j = roi.x1; j <= roi.x2; j++) {
            dm[i * width + j].hor = width;
            dm[i * width + j].ver = height;
        }
    }
    for (int i = roi.x1; i <= roi.x2; i++) {
        /// first row
        if (im[roi.y1 * width + i] == seg) {
            dm[roi.y1 * width + i].hor = 0;
            dm[roi.y1 * width + i].ver = 0;
        }
        /// last row
        if (im[roi.y2 * width + i] == seg) {
            dm[roi.y2 * width + i].hor = 0;
            dm[roi.y2 * width + i].ver = 0;
        }
    }
    for (int i = roi.y1; i <= roi.y2; i++) {
        /// first column
        if (im[i * width + roi.x1] == seg) {
            dm[i * width + roi.x1].hor = 0;
            dm[i * width + roi.x1].ver = 0;
        }
        /// last column
        if (im[i * width + roi.x2] == seg) {
            dm[i * width + roi.x2].hor = 0;
            dm[i * width + roi.x2].ver = 0;
        }
    }
    /// rest of image
    for (int i = roi.y1 + 1; i < roi.y2; i++) {
        for (int j = roi.x1 + 1; j < roi.x2; j++) {
            if ((im[i * width + j] == seg) && ((im[i * width + j - 1] != seg) || (im[i * width + j + 1] != seg) ||
                                               (im[(i - 1) * width + j] != seg) || (im[(i + 1) * width + j] != seg))) {
                dm[i * width + j].hor = 0;
                dm[i * width + j].ver = 0;
            }
        }
    }
}

/// computes normalized distance between borders of two segments from two segmentations
float borderDistance(pixtype* im1, pixtype* im2, int height, int width, pixtype seg1, pixtype seg2) {
    /// find roi
    dmroi roi;
    roi.x1 = width - 1;
    roi.y1 = height - 1;
    roi.x2 = 0;
    roi.y2 = 0;

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            if ((im1[i * width + j] == seg1) || (im2[i * width + j] == seg2)) {
                roi.x1 = min(roi.x1, j);
                roi.y1 = min(roi.y1, i);
                roi.x2 = max(roi.x2, j);
                roi.y2 = max(roi.y2, i);
            }
        }
    }

    if (!((roi.x2 - roi.x1 >= 0) && (roi.y2 - roi.y1 >= 0))) return 0;

    dmpixtype* dm1;  /// distance map
    dmpixtype* dm2;  /// distance map

    dm1 = (dmpixtype*)malloc(height * width * sizeof(dmpixtype));
#ifdef EXPORT_DISTANCE_MAP
    memset(dm1, 0, height * width * sizeof(dmpixtype));
#endif
    createBorderDistanceMapROI(im1, dm1, height, width, seg1, roi);
    computeDistanceMapROI(dm1, height, width, roi);

#ifdef EXPORT_DISTANCE_MAP
    if (0 <= seg1 && seg1 < 26) {
        char filename1[20] = "1A.raw";
        filename1[1] = (char)(seg1 + 65);
        export_distance_map(filename1, dm1, width, height);
    }
#endif

    dm2 = (dmpixtype*)malloc(height * width * sizeof(dmpixtype));
#ifdef EXPORT_DISTANCE_MAP
    memset(dm2, 0, height * width * sizeof(dmpixtype));
#endif
    createBorderDistanceMapROI(im2, dm2, height, width, seg2, roi);
    computeDistanceMapROI(dm2, height, width, roi);

#ifdef EXPORT_DISTANCE_MAP
    if (0 <= seg2 && seg2 < 26) {
        char filename2[20] = "2A.raw";
        filename2[1] = (char)(seg2 + 65);
        export_distance_map(filename2, dm2, width, height);
    }
#endif

    float sum1 = 0, sum2 = 0;
    int bord1 = 0, bord2 = 0;
    for (int i = 0; i < height * width; i++) {
        if ((dm1[i].hor == 0) && (dm1[i].ver == 0)) {  /// border pixel
            sum1 += (float)euclid(dm2[i].hor, dm2[i].ver);
            bord1++;
        }
        if ((dm2[i].hor == 0) && (dm2[i].ver == 0)) {
            sum2 += (float)euclid(dm1[i].hor, dm1[i].ver);
            bord2++;
        }
    }

    // sum/=2.0f;
    free(dm1);
    free(dm2);
    return sum1 / bord1 + sum2 / bord2;
}

void findGrouping121(pixtype* in1, pixtype* in2, int height, int width, int NoS1, int NoS2, int* bind1, int* bind2) {
    /// find sizes of intersections
    int* inter = (int*)malloc(NoS1 * NoS2 * sizeof(int));
    for (int i = 0; i < NoS1 * NoS2; i++) inter[i] = 0;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            inter[in1[j * width + i] * NoS2 + in2[j * width + i]]++;
        }
    }

    int* processed1 = (int*)malloc(NoS1 * sizeof(int));
    int* processed2 = (int*)malloc(NoS2 * sizeof(int));
    for (int i = 0; i < NoS1; i++) processed1[i] = 0;
    for (int i = 0; i < NoS2; i++) processed2[i] = 0;
    for (int i = 0; i < NoS1; i++) bind1[i] = -1;
    for (int i = 0; i < NoS2; i++) bind2[i] = -1;

    int nob = 0;  /// number of bindings
    /// find bindings
    while (1) {
        /// find max nonprocessed intersection
        int maxi = 0, maxj = 0, max = 0;
        for (int i = 1; i < NoS1; i++) {      //< ignore segment 0
            for (int j = 1; j < NoS2; j++) {  //< ignore segment 0
                if (processed1[i] || processed2[j]) continue;
                if (inter[i * NoS2 + j] > max) {
                    maxi = i;
                    maxj = j;
                    max = inter[i * NoS2 + j];
                }
            }
        }
        if (max == 0) break;
        bind1[maxi] = nob;
        bind2[maxj] = nob;
        nob++;
        processed1[maxi] = 1;
        processed2[maxj] = 1;
    }
    free(inter);
    free(processed1);
    free(processed2);
}

void implicitGrouping(int NoS1, int NoS2, int* bind1, int* bind2) {
    for (int i = 0; i < NoS1; i++) bind1[i] = i;
    for (int i = 0; i < NoS2; i++) bind2[i] = i;
}

float SD3Border(int NoS1, int NoS2, pixtype* in1, pixtype* in2, int width, int height, int* bind1, int* bind2, int nob) {
    /// merge segments

    /// stores indexes of first segment from current binding
    /// with tale element for faster searching
    int* set1 = (int*)malloc(NoS1 * sizeof(int));
    int* set2 = (int*)malloc(NoS2 * sizeof(int));

    for (int i = 0; i < NoS1; i++) set1[i] = i;
    for (int i = 0; i < NoS2; i++) set2[i] = i;

    /// fill sets
    int first = 0;
    int j = 0;
    for (int i = 0; i < nob; i++) {
        j = 0;
        bind1[NoS1] = i;  /// set tail element
        while (bind1[j] != i) j++;
        if (j < NoS1) {
            first = j;
            for (j = first; j < NoS1; j++) {
                if (bind1[j] == i) set1[j] = first;
            }
        }
        j = 0;
        bind2[NoS2] = i;
        while (bind2[j] != i) j++;
        if (j < NoS2) {
            first = j;
            for (j = first; j < NoS2; j++) {
                if (bind2[j] == i) set2[j] = first;
            }
        }
    }

    /// rewrite images - merge each binding into a single segment
    for (int i = 0; i < height * width; i++) {
        in1[i] = set1[in1[i]];
        in2[i] = set2[in2[i]];
    }
    free(set1);
    free(set2);

    /// get sizes of segments and their intersections
    // int *inter = (int*)malloc(NoS1*NoS2*sizeof(int));
    int* sizeseg1 = (int*)malloc(NoS1 * sizeof(int));
    int* sizeseg2 = (int*)malloc(NoS2 * sizeof(int));

    /// recompute sizes of segments
    for (int j = 0; j < NoS1; j++) sizeseg1[j] = 0;
    for (int j = 0; j < NoS2; j++) sizeseg2[j] = 0;

    for (int j = 0; j < height; j++) {
        for (int k = 0; k < width; k++) {
            sizeseg1[in1[j * width + k]]++;
            sizeseg2[in2[j * width + k]]++;
        }
    }

    /// for each binding compute border distance
    float bordersum = 0.0f;
    float segment_sum = 0.0f;
    for (int i = 0; i < nob; i++) {
        int j = 0;
        bind1[NoS1] = i;
        while (bind1[j] != i) j++;
        int k = 0;
        bind2[NoS2] = i;
        while (bind2[k] != i) k++;
        if ((j < NoS1) && (k < NoS2) && (sizeseg1[j] > 0) && (sizeseg2[k] > 0)) {
            segment_sum += sizeseg1[j] + sizeseg2[k];
            bordersum += (sizeseg1[j] + sizeseg2[k]) * borderDistance(in1, in2, height, width, j, k);
        }
    }

    /// extract null segments
    int sizenullseg = 0;
    for (int i = 1; i < NoS1; i++) {  //< ignore segment 0
        if (bind1[i] == -1) {
            sizenullseg += sizeseg1[i];
            segment_sum += sizeseg1[i];
        }
    }
    for (int i = 1; i < NoS2; i++) {  //< ignore segment 0
        if (bind2[i] == -1) {
            sizenullseg += sizeseg2[i];
            segment_sum += sizeseg2[i];
        }
    }
    bordersum += sizenullseg;

    bordersum = bordersum * 2.0f / (3.0f * sqrt((float)height * width) * segment_sum);

    free(sizeseg1);
    free(sizeseg2);

    return bordersum;
}

/// SD with one-many & border distance
float SD3(int NoS1, int NoS2, pixtype* in1, pixtype* in2, int width, int height) {
    /// value represents number of binding for each segment
    /// with tale element for faster searching
    int* bind1 = (int*)malloc((NoS1 + 1) * sizeof(int));
    int* bind2 = (int*)malloc((NoS2 + 1) * sizeof(int));

    /// find joint segments
    // findGroupingOptim(in1, in2, height, width, NoS1, NoS2, bind1, bind2);
    // findGrouping121(in1, in2, height, width, NoS1, NoS2, bind1, bind2);
    implicitGrouping(NoS1, NoS2, bind1, bind2);

    int nob = 0;  /// number of bindings
    for (int i = 0; i < NoS1; i++) nob = max(nob, bind1[i]);
    nob++;

    float bordersum = SD3Border(NoS1, NoS2, in1, in2, width, height, bind1, bind2, nob);

    free(bind1);
    free(bind2);
    // free(inter);
    return bordersum;
}

// renumbers segments densely from 0
void rewriteSegments(pixtype* img, int height, int width) {
    int nrew = width * height + 1;
    uint32_t* rew = (uint32_t*)malloc(nrew * sizeof(uint32_t));
    int i, j;

    for (i = 0; i < nrew; i++) rew[i] = -1;
    j = 0;
    for (i = 0; i < width * height; i++) {  /// find new numbering
        if (rew[j] != img[i]) {
            j = 0;
            while ((rew[j] != -1) && (rew[j] != img[i])) j++;
            if (rew[j] == -1) rew[j] = img[i];  ///? new segment
        }
    }
    j = 0;
    for (i = 0; i < width * height; i++) {  /// rewrite image
        /// find segment
        if (rew[j] != img[i]) {
            j = 0;
            while (rew[j] != img[i]) j++;
        }
        img[i] = j;
    }
    free(rew);
}

float callmethod(pixtype* segm1, pixtype* segm2, int width, int height) {
    // rewriteSegments(segm1, height, width);
    // rewriteSegments(segm2, height, width);

    int NoS1 = 0, NoS2 = 0;  /// number of segments
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            NoS1 = max(uint32_t(NoS1), segm1[j * width + i]);
            NoS2 = max(uint32_t(NoS2), segm2[j * width + i]);
        }
    }
    NoS1++;
    NoS2++;

    return SD3(NoS1, NoS2, segm1, segm2, width, height);
}

void processSegmentations(int width, int height, char* file1, char* file2) {
    pixtype* segm1 = (pixtype*)malloc(width * height * sizeof(pixtype));
    pixtype* segm2 = (pixtype*)malloc(width * height * sizeof(pixtype));
    pixtype* tmp = (pixtype*)malloc(width * height * sizeof(pixtype));
    FILE* fileseg;
    // FILE *ferr=_fsopen("error.txt","wb",SH_DENYWR);

    /// read 1st file
    fileseg = _fsopen(file1, "rb", SH_DENYWR);
    if (fileseg == 0) {
        cout << "Cannot open " << file1;
        return;
    }

    fread(segm1, sizeof(pixtype), width * height, fileseg);
    fclose(fileseg);

    /// read 2nd file
    fileseg = _fsopen(file2, "rb", SH_DENYWR);
    if (fileseg == 0) {
        cout << "Cannot open " << file2;
        return;
    }

    fread(segm2, sizeof(pixtype), width * height, fileseg);
    fclose(fileseg);

    memcpy(tmp, segm1, width * height * sizeof(pixtype));
    float diff = callmethod(tmp, segm2, width, height);
    cout << diff;
    return;
}

int main(int argc, char* argv[]) {
    if (argc != 5) cout << "Provide width, height, and 2 paths to files with segmentations (RAW, 16 bit).";
    int width = atoi(argv[1]);
    int height = atoi(argv[2]);
    if (width <= 0) cout << "Invalid width.";
    if (height <= 0) cout << "Invalid height.";

    processSegmentations(width, height, argv[3], argv[4]);
    return 0;
}