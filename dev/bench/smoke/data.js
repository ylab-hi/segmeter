window.BENCHMARK_DATA = {
  "lastUpdate": 1792028548191,
  "repoUrl": "https://github.com/ylab-hi/segmeter",
  "entries": {
    "Smoke (relative to bedtools)": [
      {
        "commit": {
          "author": {
            "email": "richard.schaefer@northwestern.edu",
            "name": "Richard A. Schäfer",
            "username": "riasc"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "c6fd3553435eb458eae6426d4ff8709a02b44616",
          "message": "Track continuous benchmarks and run full workloads on Slurm (#78)\n\n* Track benchmark regressions in CI and support reproducible Slurm runs\n\n* Avoid duplicate smoke benchmarks on PR branch pushes\n\n* Keep full benchmarks on Slurm and archive results manually\n\n* Upload benchmark reports without traversing private tool indexes\n\n* Require supported host Python and load its cluster module\n\n* Use repeated smoke measurements to distinguish timing noise\n\n* Run smoke benchmarks at 10K only to fit the job timeout\n\nThree repeats at 100K took about 95 minutes per source, mostly awk's\n108 s per query, so base plus head exceeded the 120-minute limit.\n\n* Record installed packages without failing on older images\n\nPublished images such as giggle-v0.14.1 bundle htslib and lack\nlibhts-dev, so dpkg-query exited 1 and aborted the cluster run before\nany benchmark started. Keep dpkg's report, including missing packages,\nas provenance instead.\n\n* Install gia and granges where Singularity users can run them\n\ncargo put the binaries under /root/.cargo/bin, and /root is 0700, so\nthe rust-tools image works under Docker (root) but not under\nSingularity as an unprivileged user: the cluster run failed with\n'Missing tool: gia'. Install to /usr/local/bin instead, and document\nthat release images up to v0.14.1 cannot run the cluster benchmark.\n\n* Back the container's /tmp with disk in Singularity runs\n\n--containall mounts a tmpfs /tmp capped by the site's session directory\nlimit (64 MiB by default), and segmeter writes each tool's output there.\nOn the Zenodo data the 100K complex 10% bin yields 77 MB from bedtools,\nso the file was truncated mid-line and scoring failed on a two-field\nline. Pass --workdir OUTPUT/tmp so /tmp lives on disk.\n\n* Document continuous benchmarking in the changelog",
          "timestamp": "2026-10-14T17:05:23-05:00",
          "tree_id": "bd7b8959aeef2c9e845388dc0232656206988cfa",
          "url": "https://github.com/ylab-hi/segmeter/commit/c6fd3553435eb458eae6426d4ff8709a02b44616"
        },
        "date": 1792028547862,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "ailist/10K/basic/3p-partial_100%",
            "value": 0.5448315248429468,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00954, 0.00934, 0.00954], \"comparison_samples\": [0.5511265164644714, 0.5334094802969731, 0.5420454545454545]}"
          },
          {
            "name": "ailist/10K/basic/5p-partial_100%",
            "value": 0.5503432494279177,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01077, 0.00943, 0.00962], \"comparison_samples\": [0.5904605263157895, 0.5394736842105263, 0.5535097813578826]}"
          },
          {
            "name": "ailist/10K/basic/contained_100%",
            "value": 0.5152838427947599,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00948, 0.00937, 0.00944], \"comparison_samples\": [0.5243362831858408, 0.5029522275899088, 0.5152838427947599]}"
          },
          {
            "name": "ailist/10K/basic/enclosed_100%",
            "value": 0.4899432697266632,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00942, 0.0095, 0.00952], \"comparison_samples\": [0.48581743166580704, 0.4945340968245705, 0.44757874941231784]}"
          },
          {
            "name": "ailist/10K/basic/left-adjacent-gap_100%",
            "value": 0.611254851228978,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00945, 0.00921, 0.0096], \"comparison_samples\": [0.6160365058670144, 0.5957309184993531, 0.6003752345215759]}"
          },
          {
            "name": "ailist/10K/basic/mid-gap1_100%",
            "value": 0.6012903225806452,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00932, 0.00923, 0.00946], \"comparison_samples\": [0.5846925972396487, 0.5962532299741602, 0.6103225806451613]}"
          },
          {
            "name": "ailist/10K/basic/mid-gap2_100%",
            "value": 0.6008997429305913,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00935, 0.00922, 0.00951], \"comparison_samples\": [0.6008997429305913, 0.6118115461181155, 0.6072796934865899]}"
          },
          {
            "name": "ailist/10K/basic/perfect-gap_100%",
            "value": 0.5913911416094823,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00947, 0.00948, 0.00958], \"comparison_samples\": [0.6093951093951093, 0.5913911416094823, 0.5928217821782178]}"
          },
          {
            "name": "ailist/10K/basic/perfect_100%",
            "value": 0.4832750873689466,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00968, 0.0096, 0.00977], \"comparison_samples\": [0.47381302006852666, 0.47928107838242634, 0.5046487603305785]}"
          },
          {
            "name": "ailist/10K/basic/right-adjacent-gap_100%",
            "value": 0.6002554278416348,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00941, 0.0092, 0.0094], \"comparison_samples\": [0.5974603174603175, 0.6013071895424836, 0.6002554278416348]}"
          },
          {
            "name": "ailist/10K/complex/mult_100bin",
            "value": 0.4897959183673469,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00593, 0.006, 0.00605], \"comparison_samples\": [0.4900826446280992, 0.4878048780487805, 0.49387755102040815]}"
          },
          {
            "name": "awk/10K/basic/3p-partial_100%",
            "value": 86.10736721873215,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.47825, 1.5808, 1.50774], \"comparison_samples\": [85.39861351819758, 90.27984009137634, 85.66704545454546]}"
          },
          {
            "name": "awk/10K/basic/5p-partial_100%",
            "value": 85.61899313501145,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.49272, 1.49956, 1.49662], \"comparison_samples\": [81.83771929824562, 85.78718535469108, 86.11162255466053]}"
          },
          {
            "name": "awk/10K/basic/contained_100%",
            "value": 81.50927947598254,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.48879, 1.49325, 1.49926], \"comparison_samples\": [82.3445796460177, 80.15297906602254, 81.83733624454149]}"
          },
          {
            "name": "awk/10K/basic/enclosed_100%",
            "value": 76.91387313047963,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.49136, 1.48605, 1.51133], \"comparison_samples\": [76.91387313047963, 77.35814679854242, 71.054536906441]}"
          },
          {
            "name": "awk/10K/basic/left-adjacent-gap_100%",
            "value": 96.23544631306598,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.49297, 1.4878, 1.4801], \"comparison_samples\": [97.32529335071708, 96.23544631306598, 92.56410256410255]}"
          },
          {
            "name": "awk/10K/basic/mid-gap1_100%",
            "value": 96.66,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.53086, 1.48711, 1.49823], \"comparison_samples\": [96.03889585947303, 96.06653746770024, 96.66]}"
          },
          {
            "name": "awk/10K/basic/mid-gap2_100%",
            "value": 95.75385604113112,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.48993, 1.48962, 1.49759], \"comparison_samples\": [95.75385604113112, 98.84671532846716, 95.63154533844188]}"
          },
          {
            "name": "awk/10K/basic/perfect-gap_100%",
            "value": 92.44853399875235,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.48195, 1.47753, 1.49638], \"comparison_samples\": [95.36357786357787, 92.17280099812851, 92.59777227722772]}"
          },
          {
            "name": "awk/10K/basic/perfect_100%",
            "value": 73.75187219171244,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.47333, 1.47725, 1.50605], \"comparison_samples\": [72.11600587371512, 73.75187219171244, 77.79183884297522]}"
          },
          {
            "name": "awk/10K/basic/right-adjacent-gap_100%",
            "value": 95.73180076628353,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.49916, 1.49806, 1.54318], \"comparison_samples\": [95.18476190476191, 97.9124183006536, 98.54278416347381]}"
          },
          {
            "name": "awk/10K/complex/mult_100bin",
            "value": 5.902857142857142,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.07231, 0.0736, 0.07165], \"comparison_samples\": [5.9760330578512395, 5.983739837398374, 5.848979591836735]}"
          },
          {
            "name": "bedmaps/10K/basic/3p-partial_100%",
            "value": 1.5231296402055967,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0269, 0.02667, 0.02656], \"comparison_samples\": [1.5540150202195264, 1.5231296402055967, 1.509090909090909]}"
          },
          {
            "name": "bedmaps/10K/basic/5p-partial_100%",
            "value": 1.5217391304347827,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02709, 0.0266, 0.0264], \"comparison_samples\": [1.4851973684210527, 1.5217391304347827, 1.518987341772152]}"
          },
          {
            "name": "bedmaps/10K/basic/contained_100%",
            "value": 1.4574235807860263,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02723, 0.02664, 0.0267], \"comparison_samples\": [1.5060840707964602, 1.429951690821256, 1.4574235807860263]}"
          },
          {
            "name": "bedmaps/10K/basic/enclosed_100%",
            "value": 1.3888602372356884,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02729, 0.02683, 0.02693], \"comparison_samples\": [1.4074265085095408, 1.3966684018740239, 1.2661024917724493]}"
          },
          {
            "name": "bedmaps/10K/basic/left-adjacent-gap_100%",
            "value": 1.4391979301423026,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02233, 0.02196, 0.02225], \"comparison_samples\": [1.455671447196871, 1.4204398447606728, 1.3914946841776108]}"
          },
          {
            "name": "bedmaps/10K/basic/mid-gap1_100%",
            "value": 1.4270967741935485,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02214, 0.02143, 0.02212], \"comparison_samples\": [1.3889585947302385, 1.3843669250645996, 1.4270967741935485]}"
          },
          {
            "name": "bedmaps/10K/basic/mid-gap2_100%",
            "value": 1.417737789203085,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02215, 0.02167, 0.02206], \"comparison_samples\": [1.423521850899743, 1.437956204379562, 1.408684546615581]}"
          },
          {
            "name": "bedmaps/10K/basic/perfect-gap_100%",
            "value": 1.394260761072988,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02238, 0.02182, 0.02235], \"comparison_samples\": [1.4401544401544402, 1.3611977542108546, 1.3830445544554453]}"
          },
          {
            "name": "bedmaps/10K/basic/perfect_100%",
            "value": 1.3564653020469297,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0276, 0.02698, 0.02717], \"comparison_samples\": [1.3509544787077827, 1.3469795307039443, 1.403409090909091]}"
          },
          {
            "name": "bedmaps/10K/basic/right-adjacent-gap_100%",
            "value": 1.4080459770114941,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02239, 0.02185, 0.02205], \"comparison_samples\": [1.4215873015873015, 1.428104575163399, 1.4080459770114941]}"
          },
          {
            "name": "bedmaps/10K/complex/mult_100bin",
            "value": 1.4685714285714284,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01799, 0.01754, 0.01824], \"comparison_samples\": [1.4867768595041322, 1.4260162601626016, 1.4889795918367346]}"
          },
          {
            "name": "bedops/10K/basic/3p-partial_100%",
            "value": 1.4431753283837807,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02527, 0.02485, 0.02586], \"comparison_samples\": [1.4598497978047373, 1.4191890348372358, 1.4693181818181817]}"
          },
          {
            "name": "bedops/10K/basic/5p-partial_100%",
            "value": 1.445080091533181,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02518, 0.02526, 0.02578], \"comparison_samples\": [1.380482456140351, 1.445080091533181, 1.4833141542002302]}"
          },
          {
            "name": "bedops/10K/basic/contained_100%",
            "value": 1.375,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.025, 0.02533, 0.02519], \"comparison_samples\": [1.3827433628318586, 1.3596349973161566, 1.375]}"
          },
          {
            "name": "bedops/10K/basic/enclosed_100%",
            "value": 1.300154718927282,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02515, 0.02521, 0.02566], \"comparison_samples\": [1.2970603403816399, 1.312337324310255, 1.2063939821344616]}"
          },
          {
            "name": "bedops/10K/basic/left-adjacent-gap_100%",
            "value": 1.4605433376455368,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02221, 0.02298, 0.02258], \"comparison_samples\": [1.4478487614080835, 1.4864165588615783, 1.41213258286429]}"
          },
          {
            "name": "bedops/10K/basic/mid-gap1_100%",
            "value": 1.4322580645161291,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0222, 0.02206, 0.02241], \"comparison_samples\": [1.3927227101631119, 1.425064599483204, 1.4458064516129032]}"
          },
          {
            "name": "bedops/10K/basic/mid-gap2_100%",
            "value": 1.4293059125964012,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02222, 0.02224, 0.02235], \"comparison_samples\": [1.4280205655526994, 1.4757796947577968, 1.4272030651340994]}"
          },
          {
            "name": "bedops/10K/basic/perfect-gap_100%",
            "value": 1.3824079850280724,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02201, 0.02216, 0.02254], \"comparison_samples\": [1.4163449163449162, 1.3824079850280724, 1.3948019801980198]}"
          },
          {
            "name": "bedops/10K/basic/perfect_100%",
            "value": 1.275586620069895,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02524, 0.02555, 0.02596], \"comparison_samples\": [1.2354380812530592, 1.275586620069895, 1.340909090909091]}"
          },
          {
            "name": "bedops/10K/basic/right-adjacent-gap_100%",
            "value": 1.434227330779055,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02193, 0.02246, 0.0227], \"comparison_samples\": [1.3923809523809525, 1.4679738562091504, 1.4495530012771394]}"
          },
          {
            "name": "bedops/10K/complex/mult_100bin",
            "value": 1.4767346938775507,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01809, 0.01806, 0.01814], \"comparison_samples\": [1.4950413223140495, 1.4682926829268292, 1.4808163265306122]}"
          },
          {
            "name": "bedtk/10K/basic/3p-partial_100%",
            "value": 0.4003426613363792,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00717, 0.00699, 0.00701], \"comparison_samples\": [0.4142114384748701, 0.3992004568817818, 0.3982954545454545]}"
          },
          {
            "name": "bedtk/10K/basic/5p-partial_100%",
            "value": 0.4101830663615561,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00699, 0.00732, 0.00717], \"comparison_samples\": [0.3832236842105263, 0.4187643020594966, 0.4125431530494822]}"
          },
          {
            "name": "bedtk/10K/basic/contained_100%",
            "value": 0.3908296943231441,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00705, 0.00721, 0.00716], \"comparison_samples\": [0.3899336283185841, 0.3870101986044015, 0.3908296943231441]}"
          },
          {
            "name": "bedtk/10K/basic/enclosed_100%",
            "value": 0.3682310469314079,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00688, 0.00714, 0.00732], \"comparison_samples\": [0.3548220732336256, 0.3716814159292035, 0.3441466854724965]}"
          },
          {
            "name": "bedtk/10K/basic/left-adjacent-gap_100%",
            "value": 0.3835705045278137,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0059, 0.00604, 0.00593], \"comparison_samples\": [0.38461538461538464, 0.39068564036222514, 0.3708567854909318]}"
          },
          {
            "name": "bedtk/10K/basic/mid-gap1_100%",
            "value": 0.3709677419354839,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00578, 0.00575, 0.00575], \"comparison_samples\": [0.3626097867001255, 0.3714470284237726, 0.3709677419354839]}"
          },
          {
            "name": "bedtk/10K/basic/mid-gap2_100%",
            "value": 0.372107969151671,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00578, 0.00591, 0.00579], \"comparison_samples\": [0.3714652956298201, 0.3921698739216988, 0.36973180076628354]}"
          },
          {
            "name": "bedtk/10K/basic/perfect-gap_100%",
            "value": 0.3655645664379289,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00586, 0.00578, 0.00596], \"comparison_samples\": [0.3770913770913771, 0.3605739238927012, 0.3688118811881188]}"
          },
          {
            "name": "bedtk/10K/basic/perfect_100%",
            "value": 0.3619570644033949,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00744, 0.00723, 0.00725], \"comparison_samples\": [0.3641703377386197, 0.36095856215676486, 0.37448347107438024]}"
          },
          {
            "name": "bedtk/10K/basic/right-adjacent-gap_100%",
            "value": 0.3735632183908046,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00582, 0.00585, 0.00588], \"comparison_samples\": [0.3695238095238095, 0.3823529411764706, 0.37547892720306514]}"
          },
          {
            "name": "bedtk/10K/complex/mult_100bin",
            "value": 0.37061224489795913,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00487, 0.00447, 0.00454], \"comparison_samples\": [0.4024793388429752, 0.36341463414634145, 0.37061224489795913]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/3p-partial_100%",
            "value": 1.034266133637921,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01831, 0.01801, 0.01811], \"comparison_samples\": [1.0577700751010977, 1.0285551113649343, 1.0289772727272728]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/5p-partial_100%",
            "value": 1.0400457665903893,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0182, 0.01781, 0.01818], \"comparison_samples\": [0.9978070175438597, 1.018878718535469, 1.046029919447641]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/contained_100%",
            "value": 1.01146288209607,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01853, 0.01776, 0.01853], \"comparison_samples\": [1.0248893805309736, 0.9533011272141707, 1.01146288209607]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/enclosed_100%",
            "value": 0.9659618359979371,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01873, 0.01952, 0.01839], \"comparison_samples\": [0.9659618359979371, 1.0161374284226965, 0.8645980253878702]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/left-adjacent-gap_100%",
            "value": 1.1190168175937905,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01742, 0.01662, 0.0173], \"comparison_samples\": [1.1355932203389831, 1.07503234152652, 1.0819262038774233]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/mid-gap1_100%",
            "value": 1.0864516129032258,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01684, 0.01656, 0.01711], \"comparison_samples\": [1.0564617314930993, 1.069767441860465, 1.1038709677419356]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/mid-gap2_100%",
            "value": 1.0970437017994858,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01707, 0.01693, 0.01718], \"comparison_samples\": [1.0970437017994858, 1.1234240212342403, 1.0970625798212006]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/perfect-gap_100%",
            "value": 1.0704928259513413,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01716, 0.0167, 0.01736], \"comparison_samples\": [1.1042471042471043, 1.041796631316282, 1.0742574257425743]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/perfect_100%",
            "value": 0.9370943584623065,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01877, 0.01809, 0.01879], \"comparison_samples\": [0.9187469407733724, 0.9031452820768846, 0.9705578512396695]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/right-adjacent-gap_100%",
            "value": 1.0881226053639845,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01721, 0.0168, 0.01704], \"comparison_samples\": [1.0926984126984127, 1.0980392156862744, 1.0881226053639845]}"
          },
          {
            "name": "bedtk_sorted/10K/complex/mult_100bin",
            "value": 0.6097959183673469,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00747, 0.00709, 0.00793], \"comparison_samples\": [0.6173553719008265, 0.5764227642276423, 0.6473469387755102]}"
          },
          {
            "name": "bedtools/10K/basic/3p-partial_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01731, 0.01751, 0.0176], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/5p-partial_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01824, 0.01748, 0.01738], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/contained_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01808, 0.01863, 0.01832], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/enclosed_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01939, 0.01921, 0.02127], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/left-adjacent-gap_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01534, 0.01546, 0.01599], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/mid-gap1_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01594, 0.01548, 0.0155], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/mid-gap2_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01556, 0.01507, 0.01566], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/perfect-gap_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01554, 0.01603, 0.01616], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/perfect_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02043, 0.02003, 0.01936], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/right-adjacent-gap_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01575, 0.0153, 0.01566], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/complex/mult_100bin",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0121, 0.0123, 0.01225], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/3p-partial_100%",
            "value": 1.2992575671045117,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0224, 0.02276, 0.02275], \"comparison_samples\": [1.294049682264587, 1.2998286693318102, 1.2926136363636362]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/5p-partial_100%",
            "value": 1.2660183066361557,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02237, 0.02213, 0.02198], \"comparison_samples\": [1.2264254385964912, 1.2660183066361557, 1.2646720368239355]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/contained_100%",
            "value": 1.237991266375546,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02269, 0.02268, 0.02257], \"comparison_samples\": [1.2549778761061947, 1.217391304347826, 1.2319868995633187]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/enclosed_100%",
            "value": 1.1712222795255287,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02305, 0.02271, 0.02229], \"comparison_samples\": [1.1887570912841672, 1.1821967725143154, 1.0479548660084625]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/left-adjacent-gap_100%",
            "value": 1.3078913324708925,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02011, 0.0204, 0.02022], \"comparison_samples\": [1.3109517601043026, 1.3195342820181113, 1.2645403377110693]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/mid-gap1_100%",
            "value": 1.2851612903225806,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02004, 0.01992, 0.01973], \"comparison_samples\": [1.2572145545796738, 1.2868217054263567, 1.2729032258064517]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/mid-gap2_100%",
            "value": 1.2769922879177378,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02032, 0.01987, 0.01982], \"comparison_samples\": [1.3059125964010285, 1.3185136031851359, 1.2656449553001277]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/perfect-gap_100%",
            "value": 1.2719900187149096,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02013, 0.02214, 0.02039], \"comparison_samples\": [1.2953667953667953, 1.3811603243917656, 1.261757425742574]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/perfect_100%",
            "value": 1.1552670993509737,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02335, 0.02314, 0.02268], \"comparison_samples\": [1.1429270680372001, 1.1552670993509737, 1.171487603305785]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/right-adjacent-gap_100%",
            "value": 1.2879948914431674,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02042, 0.02017, 0.01997], \"comparison_samples\": [1.2965079365079366, 1.3183006535947712, 1.2752234993614304]}"
          },
          {
            "name": "bedtools_sorted/10K/complex/mult_100bin",
            "value": 1.1681632653061225,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01443, 0.01431, 0.01423], \"comparison_samples\": [1.1925619834710743, 1.1634146341463414, 1.1616326530612244]}"
          },
          {
            "name": "gia/10K/basic/3p-partial_100%",
            "value": 0.9411764705882353,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01625, 0.01648, 0.0166], \"comparison_samples\": [0.9387637203928366, 0.9411764705882353, 0.9431818181818181]}"
          },
          {
            "name": "gia/10K/basic/5p-partial_100%",
            "value": 0.9485125858123571,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01729, 0.01658, 0.01653], \"comparison_samples\": [0.9479166666666667, 0.9485125858123571, 0.9510932105868815]}"
          },
          {
            "name": "gia/10K/basic/contained_100%",
            "value": 0.908296943231441,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01664, 0.01645, 0.01722], \"comparison_samples\": [0.9203539823008849, 0.8829844337090713, 0.9399563318777292]}"
          },
          {
            "name": "gia/10K/basic/enclosed_100%",
            "value": 0.8679731820526044,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01726, 0.01683, 0.01649], \"comparison_samples\": [0.890149561629706, 0.8761061946902655, 0.7752703338034791]}"
          },
          {
            "name": "gia/10K/basic/left-adjacent-gap_100%",
            "value": 1.016817593790427,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01591, 0.01519, 0.01572], \"comparison_samples\": [1.0371577574967406, 0.982535575679172, 0.9831144465290808]}"
          },
          {
            "name": "gia/10K/basic/mid-gap1_100%",
            "value": 0.9529032258064516,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01496, 0.01477, 0.01465], \"comparison_samples\": [0.9385194479297365, 0.9541343669250646, 0.9451612903225807]}"
          },
          {
            "name": "gia/10K/basic/mid-gap2_100%",
            "value": 0.9505141388174807,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01479, 0.01458, 0.01489], \"comparison_samples\": [0.9505141388174807, 0.9674850696748506, 0.9508301404853129]}"
          },
          {
            "name": "gia/10K/basic/perfect-gap_100%",
            "value": 0.9463505926388023,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01462, 0.01531, 0.01517], \"comparison_samples\": [0.9407979407979408, 0.9550842170929508, 0.9387376237623761]}"
          },
          {
            "name": "gia/10K/basic/perfect_100%",
            "value": 0.8851722416375437,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02577, 0.0169, 0.01773], \"comparison_samples\": [1.261380323054332, 0.8437343984023964, 0.9158057851239669]}"
          },
          {
            "name": "gia/10K/basic/right-adjacent-gap_100%",
            "value": 0.9572158365261814,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01516, 0.01499, 0.01488], \"comparison_samples\": [0.9625396825396825, 0.9797385620915033, 0.9501915708812261]}"
          },
          {
            "name": "gia/10K/complex/mult_100bin",
            "value": 0.9526530612244898,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01191, 0.01133, 0.01167], \"comparison_samples\": [0.984297520661157, 0.9211382113821138, 0.9526530612244898]}"
          },
          {
            "name": "giggle/10K/basic/3p-partial_100%",
            "value": 58.25299828669332,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.02972, 1.0163, 1.02001], \"comparison_samples\": [59.487001733102254, 58.0411193603655, 57.955113636363635]}"
          },
          {
            "name": "giggle/10K/basic/5p-partial_100%",
            "value": 58.88844393592677,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.01566, 1.02937, 1.04436], \"comparison_samples\": [55.68311403508772, 58.88844393592677, 60.0897583429229]}"
          },
          {
            "name": "giggle/10K/basic/contained_100%",
            "value": 57.00600436681222,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.04754, 1.01617, 1.04435], \"comparison_samples\": [57.9391592920354, 54.54482018250134, 57.00600436681222]}"
          },
          {
            "name": "giggle/10K/basic/enclosed_100%",
            "value": 52.5667870036101,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.01927, 1.01764, 1.02971], \"comparison_samples\": [52.5667870036101, 52.974492451848, 48.41137752703337]}"
          },
          {
            "name": "giggle/10K/basic/left-adjacent-gap_100%",
            "value": 3.0084087968952136,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.05042, 0.04651, 0.04645], \"comparison_samples\": [3.2868318122555413, 3.0084087968952136, 2.904940587867417]}"
          },
          {
            "name": "giggle/10K/basic/mid-gap1_100%",
            "value": 3.074193548387097,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04765, 0.04716, 0.05014], \"comparison_samples\": [2.9893350062735258, 3.046511627906977, 3.2348387096774194]}"
          },
          {
            "name": "giggle/10K/basic/mid-gap2_100%",
            "value": 3.1086118251928023,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04988, 0.04654, 0.04837], \"comparison_samples\": [3.205655526992288, 3.088254810882548, 3.0887611749680715]}"
          },
          {
            "name": "giggle/10K/basic/perfect-gap_100%",
            "value": 2.9089207735495943,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04639, 0.04713, 0.04663], \"comparison_samples\": [2.985199485199485, 2.9401122894572675, 2.8855198019801978]}"
          },
          {
            "name": "giggle/10K/basic/perfect_100%",
            "value": 51.01397903145283,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.01992, 1.03414, 1.02181], \"comparison_samples\": [49.92266275085658, 51.62955566650025, 52.77944214876034]}"
          },
          {
            "name": "giggle/10K/basic/right-adjacent-gap_100%",
            "value": 3.1411238825031926,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04919, 0.05233, 0.04802], \"comparison_samples\": [3.123174603174603, 3.420261437908497, 3.0664112388250317]}"
          },
          {
            "name": "giggle/10K/complex/mult_100bin",
            "value": 182.32408163265305,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [2.23266, 2.23418, 2.23347], \"comparison_samples\": [184.51735537190083, 181.64065040650405, 182.32408163265305]}"
          },
          {
            "name": "granges/10K/basic/3p-partial_100%",
            "value": 0.9851513420902341,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01684, 0.01766, 0.01725], \"comparison_samples\": [0.9728480647024842, 1.0085665334094802, 0.9801136363636364]}"
          },
          {
            "name": "granges/10K/basic/5p-partial_100%",
            "value": 0.9897025171624714,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0173, 0.01711, 0.0175], \"comparison_samples\": [0.9484649122807017, 0.9788329519450801, 1.0069044879171463]}"
          },
          {
            "name": "granges/10K/basic/contained_100%",
            "value": 0.9301310043668123,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01696, 0.01715, 0.01704], \"comparison_samples\": [0.9380530973451328, 0.920558239398819, 0.9301310043668123]}"
          },
          {
            "name": "granges/10K/basic/enclosed_100%",
            "value": 0.8942753996905621,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01734, 0.01732, 0.01751], \"comparison_samples\": [0.8942753996905621, 0.9016137428422695, 0.8232251998119418]}"
          },
          {
            "name": "granges/10K/basic/left-adjacent-gap_100%",
            "value": 0.6565329883570504,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01006, 0.01015, 0.01026], \"comparison_samples\": [0.6558018252933507, 0.6565329883570504, 0.6416510318949343]}"
          },
          {
            "name": "granges/10K/basic/mid-gap1_100%",
            "value": 0.6516129032258065,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00981, 0.0101, 0.01028], \"comparison_samples\": [0.6154328732747804, 0.6524547803617571, 0.6632258064516129]}"
          },
          {
            "name": "granges/10K/basic/mid-gap2_100%",
            "value": 0.6548843187660668,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00978, 0.01022, 0.01019], \"comparison_samples\": [0.62853470437018, 0.6781685467816855, 0.6507024265644955]}"
          },
          {
            "name": "granges/10K/basic/perfect-gap_100%",
            "value": 0.6238303181534623,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01, 0.01017, 0.00998], \"comparison_samples\": [0.6435006435006435, 0.6344354335620712, 0.6175742574257425]}"
          },
          {
            "name": "granges/10K/basic/perfect_100%",
            "value": 0.8691962056914627,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01741, 0.01754, 0.01709], \"comparison_samples\": [0.852178169358786, 0.8756864702945582, 0.8827479338842976]}"
          },
          {
            "name": "granges/10K/basic/right-adjacent-gap_100%",
            "value": 0.6494252873563219,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00993, 0.01017, 0.0103], \"comparison_samples\": [0.6304761904761904, 0.6647058823529413, 0.6577266922094508]}"
          },
          {
            "name": "granges/10K/complex/mult_100bin",
            "value": 0.9795918367346939,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.012, 0.01191, 0.01205], \"comparison_samples\": [0.9917355371900827, 0.9682926829268292, 0.983673469387755]}"
          },
          {
            "name": "igd/10K/basic/3p-partial_100%",
            "value": 1.2084523129640206,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02116, 0.02117, 0.02116], \"comparison_samples\": [1.222414789139226, 1.2090234151913193, 1.2022727272727274]}"
          },
          {
            "name": "igd/10K/basic/5p-partial_100%",
            "value": 1.1996567505720823,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02096, 0.02114, 0.02097], \"comparison_samples\": [1.1491228070175439, 1.2093821510297482, 1.2065592635212887]}"
          },
          {
            "name": "igd/10K/basic/contained_100%",
            "value": 1.1577510917030567,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02149, 0.02121, 0.0212], \"comparison_samples\": [1.1886061946902655, 1.1384863123993558, 1.1572052401746726]}"
          },
          {
            "name": "igd/10K/basic/enclosed_100%",
            "value": 1.1454357916451778,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02221, 0.02201, 0.0223], \"comparison_samples\": [1.1454357916451778, 1.1457574180114523, 1.0484250117536436]}"
          },
          {
            "name": "igd/10K/basic/left-adjacent-gap_100%",
            "value": 1.2509702457956016,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01894, 0.01975, 0.01934], \"comparison_samples\": [1.2346805736636244, 1.277490297542044, 1.2095059412132583]}"
          },
          {
            "name": "igd/10K/basic/mid-gap1_100%",
            "value": 1.2470967741935484,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01912, 0.01933, 0.01936], \"comparison_samples\": [1.1994981179422837, 1.2487080103359172, 1.2490322580645161]}"
          },
          {
            "name": "igd/10K/basic/mid-gap2_100%",
            "value": 1.231362467866324,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0191, 0.01924, 0.01916], \"comparison_samples\": [1.2275064267352185, 1.276708692767087, 1.223499361430396]}"
          },
          {
            "name": "igd/10K/basic/perfect-gap_100%",
            "value": 1.2270742358078603,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01959, 0.01967, 0.01972], \"comparison_samples\": [1.2606177606177607, 1.2270742358078603, 1.2202970297029703]}"
          },
          {
            "name": "igd/10K/basic/perfect_100%",
            "value": 1.1128307538691964,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02229, 0.02228, 0.0224], \"comparison_samples\": [1.091042584434655, 1.1123315027458813, 1.15702479338843]}"
          },
          {
            "name": "igd/10K/basic/right-adjacent-gap_100%",
            "value": 1.2420178799489143,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01968, 0.01945, 0.01939], \"comparison_samples\": [1.2495238095238095, 1.2712418300653594, 1.2381864623243934]}"
          },
          {
            "name": "igd/10K/complex/mult_100bin",
            "value": 0.8775510204081631,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01088, 0.01072, 0.01075], \"comparison_samples\": [0.8991735537190083, 0.8715447154471545, 0.8775510204081631]}"
          },
          {
            "name": "intervaltree/10K/basic/3p-partial_100%",
            "value": 18.497430039977154,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.32389, 0.31891, 0.33631], \"comparison_samples\": [18.711149624494514, 18.213021130782412, 19.108522727272724]}"
          },
          {
            "name": "intervaltree/10K/basic/5p-partial_100%",
            "value": 18.16533180778032,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.32043, 0.31464, 0.31753], \"comparison_samples\": [17.567434210526315, 18.0, 18.269850402761794]}"
          },
          {
            "name": "intervaltree/10K/basic/contained_100%",
            "value": 15.348799126637555,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.28046, 0.28119, 0.28198], \"comparison_samples\": [15.512168141592921, 15.093397745571657, 15.391921397379914]}"
          },
          {
            "name": "intervaltree/10K/basic/enclosed_100%",
            "value": 18.45951521402785,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.35793, 0.35651, 0.36146], \"comparison_samples\": [18.45951521402785, 18.55856324830817, 16.993888105312646]}"
          },
          {
            "name": "intervaltree/10K/basic/left-adjacent-gap_100%",
            "value": 17.34346701164295,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.27076, 0.26769, 0.26813], \"comparison_samples\": [17.65058670143416, 17.315006468305302, 16.768605378361475]}"
          },
          {
            "name": "intervaltree/10K/basic/mid-gap1_100%",
            "value": 17.33290322580645,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26922, 0.26665, 0.26866], \"comparison_samples\": [16.889585947302386, 17.225452196382427, 17.33290322580645]}"
          },
          {
            "name": "intervaltree/10K/basic/mid-gap2_100%",
            "value": 17.30591259640103,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.27024, 0.26642, 0.26928], \"comparison_samples\": [17.367609254498714, 17.67883211678832, 17.195402298850574]}"
          },
          {
            "name": "intervaltree/10K/basic/perfect-gap_100%",
            "value": 16.71490954460387,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26683, 0.26794, 0.27599], \"comparison_samples\": [17.170527670527672, 16.71490954460387, 17.07858910891089]}"
          },
          {
            "name": "intervaltree/10K/basic/perfect_100%",
            "value": 15.769845232151772,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.31587, 0.31408, 0.31824], \"comparison_samples\": [15.46108663729809, 15.680479281078384, 16.438016528925623]}"
          },
          {
            "name": "intervaltree/10K/basic/right-adjacent-gap_100%",
            "value": 17.16538952745849,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26881, 0.27178, 0.26778], \"comparison_samples\": [17.067301587301586, 17.76339869281046, 17.09961685823755]}"
          },
          {
            "name": "intervaltree/10K/complex/mult_100bin",
            "value": 32.02530612244898,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.39301, 0.39047, 0.39231], \"comparison_samples\": [32.4801652892562, 31.74552845528455, 32.02530612244898]}"
          },
          {
            "name": "tabix/10K/basic/3p-partial_100%",
            "value": 12.086236436322102,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21083, 0.21163, 0.21649], \"comparison_samples\": [12.179664933564414, 12.086236436322102, 12.30056818181818]}"
          },
          {
            "name": "tabix/10K/basic/5p-partial_100%",
            "value": 12.414187643020595,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21664, 0.21732, 0.217], \"comparison_samples\": [11.87719298245614, 12.432494279176202, 12.485615650172612]}"
          },
          {
            "name": "tabix/10K/basic/contained_100%",
            "value": 12.096069868995633,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.2216, 0.22162, 0.21687], \"comparison_samples\": [12.256637168141593, 11.895866881374127, 11.837882096069869]}"
          },
          {
            "name": "tabix/10K/basic/enclosed_100%",
            "value": 11.13305827746261,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21587, 0.21665, 0.21562], \"comparison_samples\": [11.13305827746261, 11.277980218636127, 10.137282557592854]}"
          },
          {
            "name": "tabix/10K/basic/left-adjacent-gap_100%",
            "value": 13.902975420439844,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21494, 0.2103, 0.21555], \"comparison_samples\": [14.01173402868318, 13.602846054333764, 13.48030018761726]}"
          },
          {
            "name": "tabix/10K/basic/mid-gap1_100%",
            "value": 13.643225806451612,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21147, 0.21352, 0.20873], \"comparison_samples\": [13.266624843161857, 13.793281653746769, 13.466451612903226]}"
          },
          {
            "name": "tabix/10K/basic/mid-gap2_100%",
            "value": 13.725578406169667,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21389, 0.20954, 0.21357], \"comparison_samples\": [13.746143958868895, 13.90444591904446, 13.63793103448276]}"
          },
          {
            "name": "tabix/10K/basic/perfect-gap_100%",
            "value": 13.54772301933874,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21805, 0.21503, 0.21717], \"comparison_samples\": [14.031531531531531, 13.414223331253899, 13.438737623762377]}"
          },
          {
            "name": "tabix/10K/basic/perfect_100%",
            "value": 10.876684972541188,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21796, 0.21658, 0.21786], \"comparison_samples\": [10.668624571708271, 10.812780828756864, 11.25309917355372]}"
          },
          {
            "name": "tabix/10K/basic/right-adjacent-gap_100%",
            "value": 13.80140485312899,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21755, 0.21497, 0.21613], \"comparison_samples\": [13.812698412698412, 14.05032679738562, 13.80140485312899]}"
          },
          {
            "name": "tabix/10K/complex/mult_100bin",
            "value": 0.9477551020408164,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01121, 0.01169, 0.01161], \"comparison_samples\": [0.9264462809917355, 0.9504065040650407, 0.9477551020408164]}"
          },
          {
            "name": "ucsc/10K/basic/3p-partial_100%",
            "value": 0.9337521416333523,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01611, 0.01635, 0.01639], \"comparison_samples\": [0.9306759098786829, 0.9337521416333523, 0.9312499999999998]}"
          },
          {
            "name": "ucsc/10K/basic/5p-partial_100%",
            "value": 0.9302059496567506,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01626, 0.01609, 0.01687], \"comparison_samples\": [0.8914473684210527, 0.9204805491990847, 0.9706559263521288]}"
          },
          {
            "name": "ucsc/10K/basic/contained_100%",
            "value": 0.8902838427947599,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01631, 0.01626, 0.01657], \"comparison_samples\": [0.9021017699115046, 0.8727858293075684, 0.9044759825327512]}"
          },
          {
            "name": "ucsc/10K/basic/enclosed_100%",
            "value": 0.8659102630221763,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01647, 0.01679, 0.01751], \"comparison_samples\": [0.8494069107787519, 0.8740239458615303, 0.8232251998119418]}"
          },
          {
            "name": "ucsc/10K/basic/left-adjacent-gap_100%",
            "value": 0.890038809831824,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01339, 0.01476, 0.01376], \"comparison_samples\": [0.8728813559322035, 0.9547218628719276, 0.8605378361475922]}"
          },
          {
            "name": "ucsc/10K/basic/mid-gap1_100%",
            "value": 0.8729032258064516,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01343, 0.01355, 0.01353], \"comparison_samples\": [0.842534504391468, 0.8753229974160206, 0.8729032258064516]}"
          },
          {
            "name": "ucsc/10K/basic/mid-gap2_100%",
            "value": 0.8624678663239075,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01365, 0.01339, 0.01342], \"comparison_samples\": [0.8772493573264782, 0.8885202388852025, 0.8569604086845466]}"
          },
          {
            "name": "ucsc/10K/basic/perfect-gap_100%",
            "value": 0.8496568933250157,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01362, 0.01359, 0.01456], \"comparison_samples\": [0.8764478764478765, 0.8477854023705552, 0.9009900990099009]}"
          },
          {
            "name": "ucsc/10K/basic/perfect_100%",
            "value": 0.871692461308038,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01767, 0.01746, 0.01742], \"comparison_samples\": [0.8649045521292218, 0.871692461308038, 0.8997933884297522]}"
          },
          {
            "name": "ucsc/10K/basic/right-adjacent-gap_100%",
            "value": 0.8665389527458494,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01355, 0.01357, 0.01397], \"comparison_samples\": [0.8603174603174603, 0.8869281045751635, 0.8920817369093231]}"
          },
          {
            "name": "ucsc/10K/complex/mult_100bin",
            "value": 1.062857142857143,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01332, 0.01286, 0.01302], \"comparison_samples\": [1.1008264462809918, 1.0455284552845527, 1.062857142857143]}"
          },
          {
            "name": "ailist/all/end-to-end/wall",
            "value": 0.7901196909999726,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.7888429209999686, 0.7901196909999726, 0.7927126719999933], \"comparison_samples\": [0.7888429209999686, 0.7901196909999726, 0.7927126719999933]}"
          },
          {
            "name": "awk/all/end-to-end/wall",
            "value": 15.705891379000036,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [15.647068624999974, 15.705891379000036, 15.759618380999996], \"comparison_samples\": [15.647068624999974, 15.705891379000036, 15.759618380999996]}"
          },
          {
            "name": "bedmaps/all/end-to-end/wall",
            "value": 0.9237751919999937,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.9231133160000127, 0.9277354710000054, 0.9237751919999937], \"comparison_samples\": [0.9231133160000127, 0.9277354710000054, 0.9237751919999937]}"
          },
          {
            "name": "bedops/all/end-to-end/wall",
            "value": 0.9091464230000383,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.9026600420000932, 0.9091464230000383, 0.9244977029999291], \"comparison_samples\": [0.9026600420000932, 0.9091464230000383, 0.9244977029999291]}"
          },
          {
            "name": "bedtk/all/end-to-end/wall",
            "value": 0.7215840170000547,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.7215840170000547, 0.7176602780000394, 0.7216871530000617], \"comparison_samples\": [0.7215840170000547, 0.7176602780000394, 0.7216871530000617]}"
          },
          {
            "name": "bedtk_sorted/all/end-to-end/wall",
            "value": 0.8458464020000065,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8493947849999586, 0.8458464020000065, 0.8446383460000106], \"comparison_samples\": [0.8493947849999586, 0.8458464020000065, 0.8446383460000106]}"
          },
          {
            "name": "bedtools/all/end-to-end/wall",
            "value": 0.8189149520000001,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8166485939999575, 0.8189149520000001, 0.8216461110000637], \"comparison_samples\": [0.8166485939999575, 0.8189149520000001, 0.8216461110000637]}"
          },
          {
            "name": "bedtools_sorted/all/end-to-end/wall",
            "value": 0.8850070630001028,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8935119479999685, 0.8842901649999249, 0.8850070630001028], \"comparison_samples\": [0.8935119479999685, 0.8842901649999249, 0.8850070630001028]}"
          },
          {
            "name": "gia/all/end-to-end/wall",
            "value": 0.7791475529999161,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.7889893439999014, 0.7791475529999161, 0.764065258999949], \"comparison_samples\": [0.7889893439999014, 0.7791475529999161, 0.764065258999949]}"
          },
          {
            "name": "giggle/all/end-to-end/wall",
            "value": 8.47984109999993,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [8.513253824000003, 8.437479032000056, 8.47984109999993], \"comparison_samples\": [8.513253824000003, 8.437479032000056, 8.47984109999993]}"
          },
          {
            "name": "granges/all/end-to-end/wall",
            "value": 0.7765199079999547,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8319765259999485, 0.7765199079999547, 0.7732190130000163], \"comparison_samples\": [0.8319765259999485, 0.7765199079999547, 0.7732190130000163]}"
          },
          {
            "name": "igd/all/end-to-end/wall",
            "value": 1.2583145109999805,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [1.251151909999976, 1.2583145109999805, 1.2622429490000968], \"comparison_samples\": [1.251151909999976, 1.2583145109999805, 1.2622429490000968]}"
          },
          {
            "name": "intervaltree/all/end-to-end/wall",
            "value": 3.9947727479999458,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [3.9947727479999458, 3.9730446059999167, 4.009126183000035], \"comparison_samples\": [3.9947727479999458, 3.9730446059999167, 4.009126183000035]}"
          },
          {
            "name": "segmeter/all/simulation/wall",
            "value": 1.3720877360000259,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [1.3720877360000259], \"comparison_samples\": [1.3720877360000259]}"
          },
          {
            "name": "tabix/all/end-to-end/wall",
            "value": 2.840598999000008,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [2.8492835189999823, 2.8335558809999384, 2.840598999000008], \"comparison_samples\": [2.8492835189999823, 2.8335558809999384, 2.840598999000008]}"
          },
          {
            "name": "ucsc/all/end-to-end/wall",
            "value": 0.8235726679999971,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8226756870000145, 0.8235726679999971, 0.8253490509999892], \"comparison_samples\": [0.8226756870000145, 0.8235726679999971, 0.8253490509999892]}"
          }
        ]
      }
    ]
  }
}