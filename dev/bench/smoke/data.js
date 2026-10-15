window.BENCHMARK_DATA = {
  "lastUpdate": 1792037347681,
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
      },
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
          "id": "1a5c22dd1183c3a5a2550d336c059ce01b38f806",
          "message": "Add genogrove, pygenogrove and polars-bio as benchmarked tools (#79)\n\n* Add genogrove, pygenogrove and polars-bio as benchmarked tools\n\nA fourth container on Ubuntu 24.04 holds them: genogrove needs GCC 13\nand polars-bio Python 3.11+. genogrove is index-based (idx -> .gg, isec -i);\nthe two Python tools are wrapper scripts like intervaltree.\n\n* Benchmark only the containers a source checkout has\n\nA PR that adds a container runs its base commit with the head's runner;\nthe base has neither the Dockerfile nor the tool choices, so the smoke\njob failed before measuring. The comparison now reports cases without a\nbaseline instead of aborting, and the compare step uses bash so its exit\nstatus survives the tee.\n\n* Document the genogrove, pygenogrove and polars-bio tools in the changelog\n\n* Pin polars-bio's threads and stack, keep chromosome names as text\n\npolars and DataFusion ran on every core while the other tools are\nsingle-threaded; polars, DataFusion and pyarrow floated within\npolars-bio's ranges; and scan_csv inferred Int64 for numbered\nchromosomes, failing on MT and X. The step-accounting test covers the\nthree new tools and the oracle test adds queries on a chromosome the\nreference lacks.\n\n* Report single-sample and sub-second slowdowns instead of failing\n\nThe first gated smoke run failed on bedtk's end-to-end wall (0.6 s\nagainst 0.9 s, container start-up jitter; its measured queries were\nwithin 3%) and on the simulation wall, which has one sample per side.\nA regression now needs two samples per side and, for rows in seconds,\nmore than a second of slowdown.",
          "timestamp": "2026-10-14T21:00:32-05:00",
          "tree_id": "9782a0ce986dc05c3b0542fe65d6b767adb85b28",
          "url": "https://github.com/ylab-hi/segmeter/commit/1a5c22dd1183c3a5a2550d336c059ce01b38f806"
        },
        "date": 1792030118584,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "ailist/10K/basic/3p-partial_100%",
            "value": 0.46559485530546624,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00713, 0.00807, 0.00724], \"comparison_samples\": [0.4585209003215434, 0.5015537600994406, 0.4689119170984456]}"
          },
          {
            "name": "ailist/10K/basic/5p-partial_100%",
            "value": 0.46688102893890676,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00714, 0.0086, 0.00726], \"comparison_samples\": [0.4591639871382637, 0.5155875299760192, 0.4875755540631296]}"
          },
          {
            "name": "ailist/10K/basic/contained_100%",
            "value": 0.4615384615384615,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00721, 0.00813, 0.00726], \"comparison_samples\": [0.4595283620140217, 0.47294938917975565, 0.4615384615384615]}"
          },
          {
            "name": "ailist/10K/basic/enclosed_100%",
            "value": 0.4614902943018159,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00737, 0.00835, 0.00722], \"comparison_samples\": [0.4614902943018159, 0.4755125284738041, 0.45437382001258647]}"
          },
          {
            "name": "ailist/10K/basic/left-adjacent-gap_100%",
            "value": 0.5828358208955223,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00781, 0.00824, 0.00708], \"comparison_samples\": [0.6359934853420196, 0.5823321554770319, 0.5283582089552239]}"
          },
          {
            "name": "ailist/10K/basic/mid-gap1_100%",
            "value": 0.5845181674565562,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00693, 0.00789, 0.0074], \"comparison_samples\": [0.54739336492891, 0.5659971305595409, 0.6171809841534612]}"
          },
          {
            "name": "ailist/10K/basic/mid-gap2_100%",
            "value": 0.5751734772552043,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00697, 0.00787, 0.00746], \"comparison_samples\": [0.5896785109983079, 0.546907574704656, 0.5751734772552043]}"
          },
          {
            "name": "ailist/10K/basic/perfect-gap_100%",
            "value": 0.5328947368421053,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00729, 0.00819, 0.00716], \"comparison_samples\": [0.5936482084690554, 0.5625, 0.5233918128654971]}"
          },
          {
            "name": "ailist/10K/basic/perfect_100%",
            "value": 0.4675398633257404,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00821, 0.00856, 0.00754], \"comparison_samples\": [0.4672737620944793, 0.4874715261958998, 0.45258103241296516]}"
          },
          {
            "name": "ailist/10K/basic/right-adjacent-gap_100%",
            "value": 0.5894495412844036,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00771, 0.0081, 0.00714], \"comparison_samples\": [0.6371900826446281, 0.5777460770328102, 0.5458715596330275]}"
          },
          {
            "name": "ailist/10K/complex/mult_100bin",
            "value": 0.47756729810568294,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00479, 0.00485, 0.00439], \"comparison_samples\": [0.512847965738758, 0.45158286778398515, 0.4376869391824526]}"
          },
          {
            "name": "awk/10K/basic/3p-partial_100%",
            "value": 71.93762057877814,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.11863, 1.14017, 1.06412], \"comparison_samples\": [71.93762057877814, 70.86202610316967, 68.91968911917098]}"
          },
          {
            "name": "awk/10K/basic/5p-partial_100%",
            "value": 72.65530546623795,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.13192, 1.12979, 1.08336], \"comparison_samples\": [72.79228295819937, 67.7332134292566, 72.75755540631296]}"
          },
          {
            "name": "awk/10K/basic/contained_100%",
            "value": 72.13922441195167,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.16413, 1.13475, 1.09478], \"comparison_samples\": [74.19566602931805, 66.01221640488656, 69.59821996185633]}"
          },
          {
            "name": "awk/10K/basic/enclosed_100%",
            "value": 70.74326862867876,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.16125, 1.12977, 1.09693], \"comparison_samples\": [72.71446462116467, 64.3376993166287, 69.03272498426682]}"
          },
          {
            "name": "awk/10K/basic/left-adjacent-gap_100%",
            "value": 84.2365671641791,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.16453, 1.12877, 1.07434], \"comparison_samples\": [94.83143322475571, 79.77173144876326, 80.17462686567164]}"
          },
          {
            "name": "awk/10K/basic/mid-gap1_100%",
            "value": 89.03317535545024,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.16582, 1.12716, 1.08903], \"comparison_samples\": [92.08688783570301, 80.85796269727403, 90.82819015846538]}"
          },
          {
            "name": "awk/10K/basic/mid-gap2_100%",
            "value": 85.82498072474942,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.14976, 1.11315, 1.07551], \"comparison_samples\": [97.27241962774956, 77.35580264072273, 82.92289899768696]}"
          },
          {
            "name": "awk/10K/basic/perfect-gap_100%",
            "value": 81.5672514619883,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.17729, 1.11584, 1.05762], \"comparison_samples\": [95.87052117263843, 76.63736263736263, 77.31140350877193]}"
          },
          {
            "name": "awk/10K/basic/perfect_100%",
            "value": 64.43906605922552,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.13461, 1.13155, 1.10665], \"comparison_samples\": [64.57655093910074, 64.43906605922552, 66.42557022809123]}"
          },
          {
            "name": "awk/10K/basic/right-adjacent-gap_100%",
            "value": 86.51605504587155,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.13163, 1.13233, 1.09036], \"comparison_samples\": [93.52314049586776, 80.76533523537803, 83.36085626911316]}"
          },
          {
            "name": "awk/10K/complex/mult_100bin",
            "value": 6.041874376869392,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.06077, 0.0606, 0.05806], \"comparison_samples\": [6.5064239828693795, 5.64245810055866, 5.788634097706879]}"
          },
          {
            "name": "bedmaps/10K/basic/3p-partial_100%",
            "value": 1.522829581993569,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02368, 0.0248, 0.0228], \"comparison_samples\": [1.522829581993569, 1.541330018645121, 1.4766839378238341]}"
          },
          {
            "name": "bedmaps/10K/basic/5p-partial_100%",
            "value": 1.5684887459807073,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02479, 0.02439, 0.02222], \"comparison_samples\": [1.5942122186495178, 1.462230215827338, 1.4922766957689724]}"
          },
          {
            "name": "bedmaps/10K/basic/contained_100%",
            "value": 1.523204068658614,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02328, 0.02533, 0.02396], \"comparison_samples\": [1.4837476099426385, 1.4735311227457824, 1.523204068658614]}"
          },
          {
            "name": "bedmaps/10K/basic/enclosed_100%",
            "value": 1.5798371947401375,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02523, 0.02457, 0.02543], \"comparison_samples\": [1.5798371947401375, 1.3992027334851938, 1.6003775959723097]}"
          },
          {
            "name": "bedmaps/10K/basic/left-adjacent-gap_100%",
            "value": 1.4358208955223881,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01924, 0.02111, 0.0188], \"comparison_samples\": [1.5667752442996743, 1.4918727915194347, 1.4029850746268657]}"
          },
          {
            "name": "bedmaps/10K/basic/mid-gap1_100%",
            "value": 1.5647709320695105,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01981, 0.02056, 0.01912], \"comparison_samples\": [1.5647709320695105, 1.4748923959827833, 1.5946622185154296]}"
          },
          {
            "name": "bedmaps/10K/basic/mid-gap2_100%",
            "value": 1.518118735543562,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01969, 0.02088, 0.018], \"comparison_samples\": [1.6658206429780031, 1.4510076441973592, 1.3878180416345411]}"
          },
          {
            "name": "bedmaps/10K/basic/perfect-gap_100%",
            "value": 1.5087719298245614,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01898, 0.02064, 0.02066], \"comparison_samples\": [1.5456026058631922, 1.4175824175824174, 1.510233918128655]}"
          },
          {
            "name": "bedmaps/10K/basic/perfect_100%",
            "value": 1.4441913439635536,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02449, 0.02577, 0.02536], \"comparison_samples\": [1.393853158793398, 1.4675398633257404, 1.5222088835534213]}"
          },
          {
            "name": "bedmaps/10K/basic/right-adjacent-gap_100%",
            "value": 1.4587155963302751,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01908, 0.02095, 0.01747], \"comparison_samples\": [1.5768595041322315, 1.494293865905849, 1.3356269113149848]}"
          },
          {
            "name": "bedmaps/10K/complex/mult_100bin",
            "value": 1.584247258225324,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01596, 0.01589, 0.01377], \"comparison_samples\": [1.708779443254818, 1.47951582867784, 1.3728813559322033]}"
          },
          {
            "name": "bedops/10K/basic/3p-partial_100%",
            "value": 1.4655948553054663,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02373, 0.02225, 0.02279], \"comparison_samples\": [1.5260450160771706, 1.3828464885021752, 1.4760362694300517]}"
          },
          {
            "name": "bedops/10K/basic/5p-partial_100%",
            "value": 1.4212218649517687,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02176, 0.0221, 0.02287], \"comparison_samples\": [1.39935691318328, 1.3249400479616307, 1.5359301544660846]}"
          },
          {
            "name": "bedops/10K/basic/contained_100%",
            "value": 1.399872854418309,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02128, 0.02274, 0.02202], \"comparison_samples\": [1.3562778840025496, 1.3228621291448517, 1.399872854418309]}"
          },
          {
            "name": "bedops/10K/basic/enclosed_100%",
            "value": 1.3825923606762678,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02274, 0.02208, 0.02147], \"comparison_samples\": [1.4239198497182215, 1.2574031890660593, 1.3511642542479545]}"
          },
          {
            "name": "bedops/10K/basic/left-adjacent-gap_100%",
            "value": 1.4186567164179102,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01901, 0.02063, 0.01896], \"comparison_samples\": [1.5480456026058633, 1.4579505300353357, 1.4149253731343283]}"
          },
          {
            "name": "bedops/10K/basic/mid-gap1_100%",
            "value": 1.522906793048973,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01947, 0.01907, 0.01928], \"comparison_samples\": [1.5379146919431281, 1.3680057388809184, 1.6080066722268556]}"
          },
          {
            "name": "bedops/10K/basic/mid-gap2_100%",
            "value": 1.5127216653816498,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01901, 0.01962, 0.0198], \"comparison_samples\": [1.6082910321489, 1.3634468380820013, 1.5265998457979955]}"
          },
          {
            "name": "bedops/10K/basic/perfect-gap_100%",
            "value": 1.5058479532163744,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01912, 0.02077, 0.0206], \"comparison_samples\": [1.5570032573289905, 1.426510989010989, 1.5058479532163744]}"
          },
          {
            "name": "bedops/10K/basic/perfect_100%",
            "value": 1.2562642369020502,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02286, 0.02206, 0.02159], \"comparison_samples\": [1.301081388730791, 1.2562642369020502, 1.2959183673469388]}"
          },
          {
            "name": "bedops/10K/basic/right-adjacent-gap_100%",
            "value": 1.5818042813455657,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02069, 0.02074, 0.01858], \"comparison_samples\": [1.709917355371901, 1.479315263908702, 1.4204892966360856]}"
          },
          {
            "name": "bedops/10K/complex/mult_100bin",
            "value": 1.5643070787637088,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01569, 0.01573, 0.014], \"comparison_samples\": [1.6798715203426124, 1.4646182495344509, 1.3958125623130608]}"
          },
          {
            "name": "bedtk/10K/basic/3p-partial_100%",
            "value": 0.3826366559485531,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00563, 0.00627, 0.00595], \"comparison_samples\": [0.3620578778135048, 0.3896830329397141, 0.38536269430051817]}"
          },
          {
            "name": "bedtk/10K/basic/5p-partial_100%",
            "value": 0.37491961414791,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00579, 0.00631, 0.00583], \"comparison_samples\": [0.3723472668810289, 0.37829736211031173, 0.39153794492948285]}"
          },
          {
            "name": "bedtk/10K/basic/contained_100%",
            "value": 0.4011443102352193,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00631, 0.00633, 0.00621], \"comparison_samples\": [0.4021669853409815, 0.3682373472949389, 0.3947870311506675]}"
          },
          {
            "name": "bedtk/10K/basic/enclosed_100%",
            "value": 0.3813400125234815,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00598, 0.0063, 0.00609], \"comparison_samples\": [0.3744520976831559, 0.3587699316628702, 0.3832599118942731]}"
          },
          {
            "name": "bedtk/10K/basic/left-adjacent-gap_100%",
            "value": 0.42388059701492536,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00571, 0.00568, 0.00548], \"comparison_samples\": [0.46498371335504884, 0.4014134275618375, 0.40895522388059696]}"
          },
          {
            "name": "bedtk/10K/basic/mid-gap1_100%",
            "value": 0.4265402843601896,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00583, 0.0054, 0.00525], \"comparison_samples\": [0.4605055292259084, 0.38737446197991393, 0.4378648874061718]}"
          },
          {
            "name": "bedtk/10K/basic/mid-gap2_100%",
            "value": 0.42020046260601385,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00545, 0.00546, 0.00493], \"comparison_samples\": [0.46108291032148896, 0.3794301598332175, 0.3801079414032382]}"
          },
          {
            "name": "bedtk/10K/basic/perfect-gap_100%",
            "value": 0.3720760233918129,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00497, 0.0053, 0.00509], \"comparison_samples\": [0.40472312703583063, 0.364010989010989, 0.3720760233918129]}"
          },
          {
            "name": "bedtk/10K/basic/perfect_100%",
            "value": 0.36674259681093396,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00644, 0.00674, 0.00601], \"comparison_samples\": [0.36653386454183273, 0.3838268792710707, 0.3607442977190876]}"
          },
          {
            "name": "bedtk/10K/basic/right-adjacent-gap_100%",
            "value": 0.41972477064220187,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00544, 0.00549, 0.00573], \"comparison_samples\": [0.44958677685950416, 0.391583452211127, 0.4380733944954128]}"
          },
          {
            "name": "bedtk/10K/complex/mult_100bin",
            "value": 0.3439680957128614,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00359, 0.00345, 0.00333], \"comparison_samples\": [0.38436830835117775, 0.32122905027932963, 0.3320039880358923]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/3p-partial_100%",
            "value": 0.9922829581993569,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01543, 0.01677, 0.01498], \"comparison_samples\": [0.9922829581993569, 1.0422622747047856, 0.9702072538860104]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/5p-partial_100%",
            "value": 0.9884244372990354,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01537, 0.01625, 0.0148], \"comparison_samples\": [0.9884244372990354, 0.9742206235011991, 0.9939556749496307]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/contained_100%",
            "value": 0.979020979020979,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0154, 0.01651, 0.01483], \"comparison_samples\": [0.9815168897386871, 0.9604421175101804, 0.9427844882390336]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/enclosed_100%",
            "value": 1.0619912335629305,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01696, 0.01702, 0.01472], \"comparison_samples\": [1.0619912335629305, 0.9692482915717541, 0.9263687853996223]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/left-adjacent-gap_100%",
            "value": 1.114179104477612,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01493, 0.0157, 0.01405], \"comparison_samples\": [1.215798045602606, 1.1095406360424027, 1.0485074626865671]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/mid-gap1_100%",
            "value": 1.217219589257504,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01541, 0.01649, 0.0146], \"comparison_samples\": [1.217219589257504, 1.1829268292682928, 1.2176814011676396]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/mid-gap2_100%",
            "value": 1.112567463377024,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01443, 0.017, 0.01406], \"comparison_samples\": [1.2208121827411167, 1.18137595552467, 1.0840400925212026]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/perfect-gap_100%",
            "value": 1.0862573099415205,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01486, 0.01543, 0.0138], \"comparison_samples\": [1.210097719869707, 1.0597527472527473, 1.0087719298245614]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/perfect_100%",
            "value": 0.9362186788154898,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01546, 0.01718, 0.01644], \"comparison_samples\": [0.8799089356858282, 0.978359908883827, 0.9867947178871548]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/right-adjacent-gap_100%",
            "value": 1.1414373088685017,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01493, 0.01561, 0.01468], \"comparison_samples\": [1.2338842975206612, 1.1134094151212555, 1.1223241590214068]}"
          },
          {
            "name": "bedtk_sorted/10K/complex/mult_100bin",
            "value": 0.5453639082751744,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0054, 0.0057, 0.00547], \"comparison_samples\": [0.5781584582441114, 0.5307262569832403, 0.5453639082751744]}"
          },
          {
            "name": "bedtools/10K/basic/3p-partial_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01555, 0.01609, 0.01544], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/5p-partial_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01555, 0.01668, 0.01489], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/contained_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01569, 0.01719, 0.01573], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/enclosed_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01597, 0.01756, 0.01589], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/left-adjacent-gap_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01228, 0.01415, 0.0134], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/mid-gap1_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01266, 0.01394, 0.01199], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/mid-gap2_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01182, 0.01439, 0.01297], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/perfect-gap_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01228, 0.01456, 0.01368], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/perfect_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01757, 0.01756, 0.01666], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/right-adjacent-gap_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0121, 0.01402, 0.01308], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/complex/mult_100bin",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00934, 0.01074, 0.01003], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/3p-partial_100%",
            "value": 1.2186495176848875,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01951, 0.01895, 0.01874], \"comparison_samples\": [1.254662379421222, 1.1777501553760101, 1.2137305699481864]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/5p-partial_100%",
            "value": 1.2823151125401928,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01856, 0.01994, 0.02019], \"comparison_samples\": [1.1935691318327974, 1.1954436450839327, 1.3559435862995297]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/contained_100%",
            "value": 1.232040686586141,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01938, 0.01955, 0.01772], \"comparison_samples\": [1.2351816443594648, 1.1372891215823153, 1.126509853782581]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/enclosed_100%",
            "value": 1.1809643080776453,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01796, 0.01886, 0.01899], \"comparison_samples\": [1.1246086412022542, 1.0740318906605921, 1.1950912523599748]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/left-adjacent-gap_100%",
            "value": 1.2417910447761193,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01651, 0.01775, 0.01664], \"comparison_samples\": [1.3444625407166124, 1.254416961130742, 1.2417910447761193]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/mid-gap1_100%",
            "value": 1.3799368088467614,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01806, 0.01747, 0.01699], \"comparison_samples\": [1.4265402843601895, 1.2532281205164992, 1.4170141784820685]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/mid-gap2_100%",
            "value": 1.2598303777949114,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01634, 0.01837, 0.01585], \"comparison_samples\": [1.3824027072758036, 1.2765809589993051, 1.2220508866615265]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/perfect-gap_100%",
            "value": 1.2097953216374269,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01655, 0.01774, 0.01655], \"comparison_samples\": [1.3477198697068404, 1.2184065934065933, 1.2097953216374269]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/perfect_100%",
            "value": 1.1685649202733486,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01869, 0.02052, 0.02098], \"comparison_samples\": [1.0637450199203187, 1.1685649202733486, 1.2593037214885954]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/right-adjacent-gap_100%",
            "value": 1.3333333333333335,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01744, 0.01787, 0.01647], \"comparison_samples\": [1.4413223140495868, 1.2746077032810272, 1.2591743119266054]}"
          },
          {
            "name": "bedtools_sorted/10K/complex/mult_100bin",
            "value": 1.0667996011964107,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0107, 0.01125, 0.01049], \"comparison_samples\": [1.145610278372591, 1.047486033519553, 1.0458624127617147]}"
          },
          {
            "name": "genogrove/10K/basic/3p-partial_100%",
            "value": 0.9607717041800643,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01494, 0.01542, 0.0148], \"comparison_samples\": [0.9607717041800643, 0.9583592293349906, 0.9585492227979274]}"
          },
          {
            "name": "genogrove/10K/basic/5p-partial_100%",
            "value": 0.9408360128617363,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01437, 0.01544, 0.01463], \"comparison_samples\": [0.9241157556270098, 0.9256594724220624, 0.9825386165211552]}"
          },
          {
            "name": "genogrove/10K/basic/contained_100%",
            "value": 0.9173553719008264,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01452, 0.01443, 0.01438], \"comparison_samples\": [0.9254302103250478, 0.8394415357766143, 0.9141767323585505]}"
          },
          {
            "name": "genogrove/10K/basic/enclosed_100%",
            "value": 0.9041953663118346,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01606, 0.0141, 0.01444], \"comparison_samples\": [1.0056355666875392, 0.8029612756264237, 0.9087476400251729]}"
          },
          {
            "name": "genogrove/10K/basic/left-adjacent-gap_100%",
            "value": 0.9485074626865672,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01271, 0.01297, 0.01262], \"comparison_samples\": [1.0350162866449513, 0.9166077738515902, 0.9417910447761193]}"
          },
          {
            "name": "genogrove/10K/basic/mid-gap1_100%",
            "value": 1.0410742496050553,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01318, 0.01254, 0.01358], \"comparison_samples\": [1.0410742496050553, 0.8995695839311335, 1.1326105087572977]}"
          },
          {
            "name": "genogrove/10K/basic/mid-gap2_100%",
            "value": 1.023130300693909,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01249, 0.01327, 0.0133], \"comparison_samples\": [1.0566835871404399, 0.9221681723419042, 1.0254433307632997]}"
          },
          {
            "name": "genogrove/10K/basic/perfect-gap_100%",
            "value": 0.9254385964912281,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01266, 0.01208, 0.0139], \"comparison_samples\": [1.0309446254071661, 0.8296703296703297, 1.0160818713450293]}"
          },
          {
            "name": "genogrove/10K/basic/perfect_100%",
            "value": 0.8866742596810935,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01564, 0.01557, 0.01546], \"comparison_samples\": [0.8901536710301652, 0.8866742596810935, 0.9279711884753901]}"
          },
          {
            "name": "genogrove/10K/basic/right-adjacent-gap_100%",
            "value": 0.9740061162079511,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01444, 0.01274, 0.01232], \"comparison_samples\": [1.193388429752066, 0.9087018544935807, 0.9418960244648318]}"
          },
          {
            "name": "genogrove/10K/complex/mult_100bin",
            "value": 1.3599202392821534,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01364, 0.01316, 0.01424], \"comparison_samples\": [1.4603854389721627, 1.2253258845437616, 1.4197407776669988]}"
          },
          {
            "name": "gia/10K/basic/3p-partial_100%",
            "value": 0.814790996784566,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01267, 0.01375, 0.01205], \"comparison_samples\": [0.814790996784566, 0.8545680546923555, 0.780440414507772]}"
          },
          {
            "name": "gia/10K/basic/5p-partial_100%",
            "value": 0.8662379421221865,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01373, 0.01347, 0.01342], \"comparison_samples\": [0.8829581993569131, 0.8075539568345323, 0.9012760241773001]}"
          },
          {
            "name": "gia/10K/basic/contained_100%",
            "value": 0.828353464717101,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01303, 0.01395, 0.01247], \"comparison_samples\": [0.8304652644996814, 0.8115183246073299, 0.7927527018436109]}"
          },
          {
            "name": "gia/10K/basic/enclosed_100%",
            "value": 0.8246712586098934,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01418, 0.01317, 0.01214], \"comparison_samples\": [0.8879148403256104, 0.75, 0.7640025173064819]}"
          },
          {
            "name": "gia/10K/basic/left-adjacent-gap_100%",
            "value": 0.8925373134328358,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01277, 0.01139, 0.01196], \"comparison_samples\": [1.0399022801302933, 0.8049469964664312, 0.8925373134328358]}"
          },
          {
            "name": "gia/10K/basic/mid-gap1_100%",
            "value": 0.9928909952606635,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01257, 0.01165, 0.01306], \"comparison_samples\": [0.9928909952606635, 0.8357245337159255, 1.0892410341951626]}"
          },
          {
            "name": "gia/10K/basic/mid-gap2_100%",
            "value": 0.8627602158828065,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01245, 0.01088, 0.01119], \"comparison_samples\": [1.0532994923857866, 0.7560806115357888, 0.8627602158828065]}"
          },
          {
            "name": "gia/10K/basic/perfect-gap_100%",
            "value": 0.8260233918128654,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01274, 0.01127, 0.0113], \"comparison_samples\": [1.0374592833876222, 0.7740384615384616, 0.8260233918128654]}"
          },
          {
            "name": "gia/10K/basic/perfect_100%",
            "value": 0.7904328018223234,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0256, 0.01388, 0.01385], \"comparison_samples\": [1.4570290267501425, 0.7904328018223234, 0.8313325330132052]}"
          },
          {
            "name": "gia/10K/basic/right-adjacent-gap_100%",
            "value": 0.9403669724770642,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01315, 0.01188, 0.0123], \"comparison_samples\": [1.0867768595041323, 0.8473609129814551, 0.9403669724770642]}"
          },
          {
            "name": "gia/10K/complex/mult_100bin",
            "value": 0.8384845463609173,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00841, 0.00752, 0.00868], \"comparison_samples\": [0.9004282655246254, 0.7001862197392924, 0.8654037886340976]}"
          },
          {
            "name": "giggle/10K/basic/3p-partial_100%",
            "value": 58.857234726688105,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.91523, 0.91866, 0.87318], \"comparison_samples\": [58.857234726688105, 57.09509011808577, 56.55310880829015]}"
          },
          {
            "name": "giggle/10K/basic/5p-partial_100%",
            "value": 56.615434083601286,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.88037, 0.91077, 0.87548], \"comparison_samples\": [56.615434083601286, 54.602517985611506, 58.79650772330423]}"
          },
          {
            "name": "giggle/10K/basic/contained_100%",
            "value": 56.95232040686586,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.91777, 0.89586, 0.86413], \"comparison_samples\": [58.49394518801785, 52.11518324607329, 54.93515575333757]}"
          },
          {
            "name": "giggle/10K/basic/enclosed_100%",
            "value": 55.80651221039448,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.922, 0.89123, 0.85788], \"comparison_samples\": [57.73324984345648, 50.75341685649203, 53.98867212083071]}"
          },
          {
            "name": "giggle/10K/basic/left-adjacent-gap_100%",
            "value": 3.2179104477611937,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04485, 0.04312, 0.0417], \"comparison_samples\": [3.65228013029316, 3.047349823321555, 3.1119402985074625]}"
          },
          {
            "name": "giggle/10K/basic/mid-gap1_100%",
            "value": 3.327804107424961,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04398, 0.04205, 0.04213], \"comparison_samples\": [3.4739336492890995, 3.0164992826398853, 3.5137614678899083]}"
          },
          {
            "name": "giggle/10K/basic/mid-gap2_100%",
            "value": 3.3392444101773324,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04531, 0.04274, 0.04331], \"comparison_samples\": [3.8333333333333335, 2.9701181375955525, 3.3392444101773324]}"
          },
          {
            "name": "giggle/10K/basic/perfect-gap_100%",
            "value": 3.2945906432748537,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04646, 0.04507, 0.04474], \"comparison_samples\": [3.7833876221498373, 3.095467032967033, 3.27046783625731]}"
          },
          {
            "name": "giggle/10K/basic/perfect_100%",
            "value": 50.43166287015946,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.88558, 0.92769, 0.8662], \"comparison_samples\": [50.40295959021059, 52.82972665148064, 51.992797118847534]}"
          },
          {
            "name": "giggle/10K/basic/right-adjacent-gap_100%",
            "value": 3.3692660550458715,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04632, 0.04407, 0.04083], \"comparison_samples\": [3.8280991735537193, 3.1433666191155494, 3.121559633027523]}"
          },
          {
            "name": "giggle/10K/complex/mult_100bin",
            "value": 188.9531405782652,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.94894, 1.8952, 1.85855], \"comparison_samples\": [208.66595289079228, 176.46182495344507, 185.29910269192422]}"
          },
          {
            "name": "granges/10K/basic/3p-partial_100%",
            "value": 0.8938906752411575,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0139, 0.01421, 0.01258], \"comparison_samples\": [0.8938906752411575, 0.8831572405220635, 0.8147668393782382]}"
          },
          {
            "name": "granges/10K/basic/5p-partial_100%",
            "value": 0.8482315112540193,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01319, 0.014, 0.01246], \"comparison_samples\": [0.8482315112540193, 0.8393285371702638, 0.8368032236400269]}"
          },
          {
            "name": "granges/10K/basic/contained_100%",
            "value": 0.8645899554990463,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01376, 0.0136, 0.01293], \"comparison_samples\": [0.8769917144678139, 0.7911576497963931, 0.8219961856325493]}"
          },
          {
            "name": "granges/10K/basic/enclosed_100%",
            "value": 0.7946149029430181,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01303, 0.01269, 0.01237], \"comparison_samples\": [0.8159048215403881, 0.7226651480637813, 0.7784770295783512]}"
          },
          {
            "name": "granges/10K/basic/left-adjacent-gap_100%",
            "value": 0.5955223880597014,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00841, 0.00798, 0.00785], \"comparison_samples\": [0.6848534201954398, 0.5639575971731449, 0.585820895522388]}"
          },
          {
            "name": "granges/10K/basic/mid-gap1_100%",
            "value": 0.631911532385466,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00819, 0.008, 0.00778], \"comparison_samples\": [0.6469194312796208, 0.5738880918220948, 0.6488740617180984]}"
          },
          {
            "name": "granges/10K/basic/mid-gap2_100%",
            "value": 0.5859676175790285,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0076, 0.00828, 0.00723], \"comparison_samples\": [0.6429780033840947, 0.5753995830437804, 0.5574402467232074]}"
          },
          {
            "name": "granges/10K/basic/perfect-gap_100%",
            "value": 0.5497076023391813,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00852, 0.00752, 0.00743], \"comparison_samples\": [0.6938110749185668, 0.5164835164835164, 0.5431286549707602]}"
          },
          {
            "name": "granges/10K/basic/perfect_100%",
            "value": 0.7761958997722096,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01378, 0.01363, 0.01314], \"comparison_samples\": [0.7842914058053501, 0.7761958997722096, 0.7887154861944777]}"
          },
          {
            "name": "granges/10K/basic/right-adjacent-gap_100%",
            "value": 0.6269113149847095,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00832, 0.00768, 0.0082], \"comparison_samples\": [0.687603305785124, 0.5477888730385164, 0.6269113149847095]}"
          },
          {
            "name": "granges/10K/complex/mult_100bin",
            "value": 0.8733798604187437,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00935, 0.00876, 0.0084], \"comparison_samples\": [1.0010706638115634, 0.8156424581005587, 0.8374875373878364]}"
          },
          {
            "name": "igd/10K/basic/3p-partial_100%",
            "value": 0.8340836012861736,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01293, 0.01297, 0.01348], \"comparison_samples\": [0.8315112540192927, 0.8060907395898074, 0.8730569948186528]}"
          },
          {
            "name": "igd/10K/basic/5p-partial_100%",
            "value": 0.8135048231511254,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01241, 0.01277, 0.01265], \"comparison_samples\": [0.7980707395498392, 0.7655875299760192, 0.8495634654130289]}"
          },
          {
            "name": "igd/10K/basic/contained_100%",
            "value": 0.7552447552447552,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01175, 0.01344, 0.01188], \"comparison_samples\": [0.7488846398980242, 0.7818499127399651, 0.7552447552447552]}"
          },
          {
            "name": "igd/10K/basic/enclosed_100%",
            "value": 0.7520350657482779,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01201, 0.0132, 0.012], \"comparison_samples\": [0.7520350657482779, 0.7517084282460137, 0.7551919446192573]}"
          },
          {
            "name": "igd/10K/basic/left-adjacent-gap_100%",
            "value": 0.8037313432835821,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01053, 0.01148, 0.01077], \"comparison_samples\": [0.8574918566775245, 0.8113074204946997, 0.8037313432835821]}"
          },
          {
            "name": "igd/10K/basic/mid-gap1_100%",
            "value": 0.8759873617693522,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01038, 0.01109, 0.0119], \"comparison_samples\": [0.8199052132701422, 0.7955523672883787, 0.9924937447873228]}"
          },
          {
            "name": "igd/10K/basic/mid-gap2_100%",
            "value": 0.8218966846569005,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01041, 0.01154, 0.01066], \"comparison_samples\": [0.8807106598984772, 0.8019457956914524, 0.8218966846569005]}"
          },
          {
            "name": "igd/10K/basic/perfect-gap_100%",
            "value": 0.807748538011696,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01105, 0.01156, 0.01087], \"comparison_samples\": [0.8998371335504887, 0.793956043956044, 0.7945906432748538]}"
          },
          {
            "name": "igd/10K/basic/perfect_100%",
            "value": 0.6867881548974943,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01206, 0.01351, 0.01205], \"comparison_samples\": [0.6863972680705749, 0.7693621867881549, 0.7232893157262905]}"
          },
          {
            "name": "igd/10K/basic/right-adjacent-gap_100%",
            "value": 0.7889908256880733,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01019, 0.01176, 0.01032], \"comparison_samples\": [0.8421487603305785, 0.8388017118402282, 0.7889908256880733]}"
          },
          {
            "name": "igd/10K/complex/mult_100bin",
            "value": 0.7048853439680957,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00688, 0.00874, 0.00707], \"comparison_samples\": [0.7366167023554604, 0.813780260707635, 0.7048853439680957]}"
          },
          {
            "name": "intervaltree/10K/basic/3p-partial_100%",
            "value": 18.194855305466238,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.2835, 0.28293, 0.27242], \"comparison_samples\": [18.231511254019292, 17.584213797389683, 17.643782383419687]}"
          },
          {
            "name": "intervaltree/10K/basic/5p-partial_100%",
            "value": 18.107395498392282,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.28558, 0.28157, 0.25755], \"comparison_samples\": [18.365273311897106, 16.880695443645084, 17.296843519140364]}"
          },
          {
            "name": "intervaltree/10K/basic/contained_100%",
            "value": 15.49713922441195,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.24377, 0.24496, 0.23525], \"comparison_samples\": [15.536647546207776, 14.250145433391507, 14.955499046408136]}"
          },
          {
            "name": "intervaltree/10K/basic/enclosed_100%",
            "value": 19.608641202254223,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.31777, 0.31315, 0.29714], \"comparison_samples\": [19.8979336255479, 17.833143507972665, 18.699811202013844]}"
          },
          {
            "name": "intervaltree/10K/basic/left-adjacent-gap_100%",
            "value": 17.301492537313433,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.23659, 0.23184, 0.22817], \"comparison_samples\": [19.26628664495114, 16.384452296819788, 17.027611940298506]}"
          },
          {
            "name": "intervaltree/10K/basic/mid-gap1_100%",
            "value": 18.419431279620856,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.2357, 0.23319, 0.21891], \"comparison_samples\": [18.617693522906794, 16.728120516499285, 18.257714762301916]}"
          },
          {
            "name": "intervaltree/10K/basic/mid-gap2_100%",
            "value": 17.855050115651505,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.23158, 0.23847, 0.22509], \"comparison_samples\": [19.592216582064296, 16.571924947880472, 17.35466461063994]}"
          },
          {
            "name": "intervaltree/10K/basic/perfect-gap_100%",
            "value": 17.111842105263158,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.23482, 0.23409, 0.21813], \"comparison_samples\": [19.122149837133552, 16.07760989010989, 15.945175438596491]}"
          },
          {
            "name": "intervaltree/10K/basic/perfect_100%",
            "value": 15.952733485193622,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.28084, 0.28013, 0.26412], \"comparison_samples\": [15.98406374501992, 15.952733485193622, 15.853541416566626]}"
          },
          {
            "name": "intervaltree/10K/basic/right-adjacent-gap_100%",
            "value": 17.43807339449541,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.22961, 0.22809, 0.22549], \"comparison_samples\": [18.97603305785124, 16.268901569186877, 17.239296636085626]}"
          },
          {
            "name": "intervaltree/10K/complex/mult_100bin",
            "value": 32.468594217347956,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.32566, 0.33234, 0.31743], \"comparison_samples\": [34.86723768736617, 30.944134078212294, 31.64805583250249]}"
          },
          {
            "name": "polars_bio/10K/basic/3p-partial_100%",
            "value": 15.752411575562702,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26028, 0.24465, 0.24495], \"comparison_samples\": [16.738263665594857, 15.205096333126166, 15.864637305699482]}"
          },
          {
            "name": "polars_bio/10K/basic/5p-partial_100%",
            "value": 16.400000000000002,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.25869, 0.25502, 0.24819], \"comparison_samples\": [16.636012861736333, 15.28896882494005, 16.668233713901948]}"
          },
          {
            "name": "polars_bio/10K/basic/contained_100%",
            "value": 15.72027972027972,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26068, 0.24728, 0.24092], \"comparison_samples\": [16.61440407903123, 14.385107620709714, 15.315956770502224]}"
          },
          {
            "name": "polars_bio/10K/basic/enclosed_100%",
            "value": 15.538509705698184,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.25978, 0.24797, 0.24815], \"comparison_samples\": [16.26675015654352, 14.121298405466971, 15.616740088105725]}"
          },
          {
            "name": "polars_bio/10K/basic/left-adjacent-gap_100%",
            "value": 19.679850746268656,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26371, 0.27157, 0.24861], \"comparison_samples\": [21.474755700325733, 19.192226148409894, 18.552985074626864]}"
          },
          {
            "name": "polars_bio/10K/basic/mid-gap1_100%",
            "value": 19.71958925750395,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26078, 0.23983, 0.24965], \"comparison_samples\": [20.59873617693523, 17.20444763271162, 20.821517931609673]}"
          },
          {
            "name": "polars_bio/10K/basic/mid-gap2_100%",
            "value": 19.164225134926753,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.27021, 0.24834, 0.24856], \"comparison_samples\": [22.860406091370557, 17.257817929117444, 19.164225134926753]}"
          },
          {
            "name": "polars_bio/10K/basic/perfect-gap_100%",
            "value": 18.937865497076025,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.25907, 0.25939, 0.24561], \"comparison_samples\": [21.096905537459286, 17.815247252747252, 17.95394736842105]}"
          },
          {
            "name": "polars_bio/10K/basic/perfect_100%",
            "value": 14.344533029612757,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.30847, 0.24344, 0.25189], \"comparison_samples\": [17.556630620375643, 13.863325740318906, 15.119447779111644]}"
          },
          {
            "name": "polars_bio/10K/basic/right-adjacent-gap_100%",
            "value": 19.448012232415902,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.2681, 0.25438, 0.24868], \"comparison_samples\": [22.15702479338843, 18.144079885877318, 19.012232415902144]}"
          },
          {
            "name": "polars_bio/10K/complex/mult_100bin",
            "value": 24.694915254237287,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.25984, 0.24769, 0.24299], \"comparison_samples\": [27.82012847965739, 23.06238361266294, 24.22632103688933]}"
          },
          {
            "name": "pygenogrove/10K/basic/3p-partial_100%",
            "value": 9.5016077170418,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.148, 0.14775, 0.14666], \"comparison_samples\": [9.517684887459806, 9.18272218769422, 9.498704663212436]}"
          },
          {
            "name": "pygenogrove/10K/basic/5p-partial_100%",
            "value": 9.479742765273313,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.15088, 0.14419, 0.14741], \"comparison_samples\": [9.70289389067524, 8.644484412470025, 9.899932840832774]}"
          },
          {
            "name": "pygenogrove/10K/basic/contained_100%",
            "value": 9.33375715193897,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.14707, 0.14682, 0.14179], \"comparison_samples\": [9.373486297004462, 8.541012216404887, 9.013986013986013]}"
          },
          {
            "name": "pygenogrove/10K/basic/enclosed_100%",
            "value": 9.103318722604884,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.14538, 0.14835, 0.14001], \"comparison_samples\": [9.103318722604884, 8.448177676537586, 8.811202013845184]}"
          },
          {
            "name": "pygenogrove/10K/basic/left-adjacent-gap_100%",
            "value": 9.044029850746268,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.12835, 0.12095, 0.12119], \"comparison_samples\": [10.451954397394136, 8.547703180212014, 9.044029850746268]}"
          },
          {
            "name": "pygenogrove/10K/basic/mid-gap1_100%",
            "value": 9.677725118483414,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.1263, 0.12252, 0.11983], \"comparison_samples\": [9.976303317535546, 8.78909612625538, 9.994161801501251]}"
          },
          {
            "name": "pygenogrove/10K/basic/mid-gap2_100%",
            "value": 9.723207401696222,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.12611, 0.12525, 0.12691], \"comparison_samples\": [10.669204737732656, 8.70396108408617, 9.784888203546645]}"
          },
          {
            "name": "pygenogrove/10K/basic/perfect-gap_100%",
            "value": 9.069444444444445,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.12486, 0.12407, 0.11871], \"comparison_samples\": [10.167752442996743, 8.521291208791208, 8.677631578947368]}"
          },
          {
            "name": "pygenogrove/10K/basic/perfect_100%",
            "value": 8.47608200455581,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.17784, 0.14884, 0.14576], \"comparison_samples\": [10.121798520204896, 8.47608200455581, 8.749099639855942]}"
          },
          {
            "name": "pygenogrove/10K/basic/right-adjacent-gap_100%",
            "value": 9.214831804281346,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.12471, 0.12031, 0.12053], \"comparison_samples\": [10.306611570247934, 8.581312410841655, 9.214831804281346]}"
          },
          {
            "name": "pygenogrove/10K/complex/mult_100bin",
            "value": 9.559322033898304,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.09588, 0.09642, 0.09489], \"comparison_samples\": [10.265524625267668, 8.977653631284918, 9.46061814556331]}"
          },
          {
            "name": "tabix/10K/basic/3p-partial_100%",
            "value": 11.187781350482316,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.17397, 0.18065, 0.17187], \"comparison_samples\": [11.187781350482316, 11.22747047855811, 11.131476683937823]}"
          },
          {
            "name": "tabix/10K/basic/5p-partial_100%",
            "value": 11.046945337620578,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.17178, 0.18052, 0.16602], \"comparison_samples\": [11.046945337620578, 10.82254196642686, 11.149764942914707]}"
          },
          {
            "name": "tabix/10K/basic/contained_100%",
            "value": 10.88239033693579,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.16873, 0.18021, 0.17118], \"comparison_samples\": [10.753983428935628, 10.483420593368237, 10.88239033693579]}"
          },
          {
            "name": "tabix/10K/basic/enclosed_100%",
            "value": 10.747651847213524,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.1648, 0.18231, 0.17164], \"comparison_samples\": [10.319348778960551, 10.382118451025057, 10.801762114537443]}"
          },
          {
            "name": "tabix/10K/basic/left-adjacent-gap_100%",
            "value": 12.302238805970148,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.1631, 0.17622, 0.16485], \"comparison_samples\": [13.281758957654723, 12.453710247349823, 12.302238805970148]}"
          },
          {
            "name": "tabix/10K/basic/mid-gap1_100%",
            "value": 13.10821484992101,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.16595, 0.17533, 0.16252], \"comparison_samples\": [13.10821484992101, 12.577474892395985, 13.55462885738115]}"
          },
          {
            "name": "tabix/10K/basic/mid-gap2_100%",
            "value": 13.026985350809559,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.1667, 0.17013, 0.16896], \"comparison_samples\": [14.103214890016918, 11.822793606671299, 13.026985350809559]}"
          },
          {
            "name": "tabix/10K/basic/perfect-gap_100%",
            "value": 12.404970760233917,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.1697, 0.1738, 0.16694], \"comparison_samples\": [13.819218241042345, 11.936813186813188, 12.203216374269006]}"
          },
          {
            "name": "tabix/10K/basic/perfect_100%",
            "value": 10.280751708428246,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.17203, 0.18053, 0.19034], \"comparison_samples\": [9.791121229368242, 10.280751708428246, 11.424969987995198]}"
          },
          {
            "name": "tabix/10K/basic/right-adjacent-gap_100%",
            "value": 12.873088685015292,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.16838, 0.17433, 0.16286], \"comparison_samples\": [13.915702479338844, 12.434379457917263, 12.451070336391439]}"
          },
          {
            "name": "tabix/10K/complex/mult_100bin",
            "value": 0.9381854436689929,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00941, 0.00892, 0.00984], \"comparison_samples\": [1.0074946466809422, 0.830540037243948, 0.9810568295114656]}"
          },
          {
            "name": "ucsc/10K/basic/3p-partial_100%",
            "value": 0.8585209003215435,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01297, 0.01403, 0.01335], \"comparison_samples\": [0.8340836012861736, 0.8719701678060908, 0.8646373056994818]}"
          },
          {
            "name": "ucsc/10K/basic/5p-partial_100%",
            "value": 0.9054662379421222,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01293, 0.01492, 0.01408], \"comparison_samples\": [0.8315112540192927, 0.8944844124700239, 0.9456010745466756]}"
          },
          {
            "name": "ucsc/10K/basic/contained_100%",
            "value": 0.8226319135410044,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01283, 0.01432, 0.01294], \"comparison_samples\": [0.8177182919056724, 0.8330424665503199, 0.8226319135410044]}"
          },
          {
            "name": "ucsc/10K/basic/enclosed_100%",
            "value": 0.8716343143393863,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01479, 0.01392, 0.01288], \"comparison_samples\": [0.9261114589855979, 0.7927107061503417, 0.8105726872246696]}"
          },
          {
            "name": "ucsc/10K/basic/left-adjacent-gap_100%",
            "value": 0.8388059701492537,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01102, 0.01169, 0.01124], \"comparison_samples\": [0.8973941368078177, 0.8261484098939931, 0.8388059701492537]}"
          },
          {
            "name": "ucsc/10K/basic/mid-gap1_100%",
            "value": 0.9060031595576619,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01124, 0.01177, 0.01147], \"comparison_samples\": [0.8878357030015798, 0.8443328550932568, 0.9566305254378648]}"
          },
          {
            "name": "ucsc/10K/basic/mid-gap2_100%",
            "value": 0.8997686969930608,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01167, 0.01189, 0.01016], \"comparison_samples\": [0.9873096446700507, 0.8262682418346073, 0.7833461835003855]}"
          },
          {
            "name": "ucsc/10K/basic/perfect-gap_100%",
            "value": 0.8238304093567252,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01096, 0.01175, 0.01127], \"comparison_samples\": [0.8925081433224755, 0.8070054945054945, 0.8238304093567252]}"
          },
          {
            "name": "ucsc/10K/basic/perfect_100%",
            "value": 0.8325740318906606,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01395, 0.01504, 0.01462], \"comparison_samples\": [0.7939669891861127, 0.856492027334852, 0.8775510204081631]}"
          },
          {
            "name": "ucsc/10K/basic/right-adjacent-gap_100%",
            "value": 0.8814984709480123,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01097, 0.0123, 0.01153], \"comparison_samples\": [0.9066115702479339, 0.877318116975749, 0.8814984709480123]}"
          },
          {
            "name": "ucsc/10K/complex/mult_100bin",
            "value": 1.11864406779661,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01155, 0.01122, 0.00979], \"comparison_samples\": [1.2366167023554604, 1.0446927374301678, 0.9760717846460617]}"
          },
          {
            "name": "ailist/all/end-to-end/wall",
            "value": 0.6810173520000262,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.6810173520000262, 0.7406484570000202, 0.6612883110000212], \"comparison_samples\": [0.6810173520000262, 0.7406484570000202, 0.6612883110000212]}"
          },
          {
            "name": "awk/all/end-to-end/wall",
            "value": 12.007326403999969,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [12.186572578999971, 12.007326403999969, 11.48994053299998], \"comparison_samples\": [12.186572578999971, 12.007326403999969, 11.48994053299998]}"
          },
          {
            "name": "bedmaps/all/end-to-end/wall",
            "value": 0.8359743629999912,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8213033060000043, 0.8670607729999915, 0.8359743629999912], \"comparison_samples\": [0.8213033060000043, 0.8670607729999915, 0.8359743629999912]}"
          },
          {
            "name": "bedops/all/end-to-end/wall",
            "value": 0.800927999999999,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.800927999999999, 0.8112026480000054, 0.7983238420000021], \"comparison_samples\": [0.800927999999999, 0.8112026480000054, 0.7983238420000021]}"
          },
          {
            "name": "bedtk/all/end-to-end/wall",
            "value": 0.6408702090000133,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.6408702090000133, 0.6754400049999845, 0.6312642269999742], \"comparison_samples\": [0.6408702090000133, 0.6754400049999845, 0.6312642269999742]}"
          },
          {
            "name": "bedtk_sorted/all/end-to-end/wall",
            "value": 0.7506150439999715,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.7506150439999715, 0.7977345840000112, 0.7256110430000149], \"comparison_samples\": [0.7506150439999715, 0.7977345840000112, 0.7256110430000149]}"
          },
          {
            "name": "bedtools/all/end-to-end/wall",
            "value": 0.7295758660000047,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.7295758660000047, 0.7647151749999921, 0.7168187259999854], \"comparison_samples\": [0.7295758660000047, 0.7647151749999921, 0.7168187259999854]}"
          },
          {
            "name": "bedtools_sorted/all/end-to-end/wall",
            "value": 0.772074210000028,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.772074210000028, 0.8017611890000467, 0.7586695710000413], \"comparison_samples\": [0.772074210000028, 0.8017611890000467, 0.7586695710000413]}"
          },
          {
            "name": "genogrove/all/end-to-end/wall",
            "value": 0.7120803390000106,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.7823117099999877, 0.7079415649999987, 0.7120803390000106], \"comparison_samples\": [0.7823117099999877, 0.7079415649999987, 0.7120803390000106]}"
          },
          {
            "name": "gia/all/end-to-end/wall",
            "value": 0.6537737839999522,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.741430078999997, 0.6537737839999522, 0.6455416560000344], \"comparison_samples\": [0.741430078999997, 0.6537737839999522, 0.6455416560000344]}"
          },
          {
            "name": "giggle/all/end-to-end/wall",
            "value": 7.414087359999996,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [7.50087166000003, 7.414087359999996, 7.214957400000003], \"comparison_samples\": [7.50087166000003, 7.414087359999996, 7.214957400000003]}"
          },
          {
            "name": "granges/all/end-to-end/wall",
            "value": 0.6628571459999648,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.9495198569999843, 0.6628571459999648, 0.6344000500000107], \"comparison_samples\": [0.9495198569999843, 0.6628571459999648, 0.6344000500000107]}"
          },
          {
            "name": "igd/all/end-to-end/wall",
            "value": 0.9312901399999873,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8909411629999795, 0.9795178300000202, 0.9312901399999873], \"comparison_samples\": [0.8909411629999795, 0.9795178300000202, 0.9312901399999873]}"
          },
          {
            "name": "intervaltree/all/end-to-end/wall",
            "value": 3.5159967120000033,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [3.5236517260000255, 3.5159967120000033, 3.3574389969999743], \"comparison_samples\": [3.5236517260000255, 3.5159967120000033, 3.3574389969999743]}"
          },
          {
            "name": "polars_bio/all/end-to-end/wall",
            "value": 3.3197303839999677,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [3.515611734999993, 3.3197303839999677, 3.2971067849999827], \"comparison_samples\": [3.515611734999993, 3.3197303839999677, 3.2971067849999827]}"
          },
          {
            "name": "pygenogrove/all/end-to-end/wall",
            "value": 2.000030038000034,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [2.063110380000012, 2.000030038000034, 1.9727177660000166], \"comparison_samples\": [2.063110380000012, 2.000030038000034, 1.9727177660000166]}"
          },
          {
            "name": "segmeter/all/simulation/wall",
            "value": 1.5039054549999946,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [1.5039054549999946], \"comparison_samples\": [1.5039054549999946]}"
          },
          {
            "name": "tabix/all/end-to-end/wall",
            "value": 2.3423034569999572,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [2.312549826999998, 2.436049910999998, 2.3423034569999572], \"comparison_samples\": [2.312549826999998, 2.436049910999998, 2.3423034569999572]}"
          },
          {
            "name": "ucsc/all/end-to-end/wall",
            "value": 0.7307997839999985,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.7307997839999985, 0.7469671609999864, 0.7104244170000129], \"comparison_samples\": [0.7307997839999985, 0.7469671609999864, 0.7104244170000129]}"
          }
        ]
      },
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
          "id": "107749720ab1d40e573b5e15100d3b0df31dd141",
          "message": "Release v0.15.0 (#80)",
          "timestamp": "2026-10-14T21:34:20-05:00",
          "tree_id": "b60df605bc42e122469f2583cfc1c4c6318719f5",
          "url": "https://github.com/ylab-hi/segmeter/commit/107749720ab1d40e573b5e15100d3b0df31dd141"
        },
        "date": 1792031979757,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "ailist/10K/basic/3p-partial_100%",
            "value": 0.5610663360198389,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00877, 0.00988, 0.00905], \"comparison_samples\": [0.5360635696821516, 0.6329276105060858, 0.5610663360198389]}"
          },
          {
            "name": "ailist/10K/basic/5p-partial_100%",
            "value": 0.5897273303741281,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0093, 0.00961, 0.00898], \"comparison_samples\": [0.5897273303741281, 0.5957842529448234, 0.5701587301587302]}"
          },
          {
            "name": "ailist/10K/basic/contained_100%",
            "value": 0.4948688711516534,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00867, 0.00868, 0.00932], \"comparison_samples\": [0.528981086028066, 0.4903954802259887, 0.5313568985176739]}"
          },
          {
            "name": "ailist/10K/basic/enclosed_100%",
            "value": 0.5008431703204047,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00897, 0.0087, 0.00891], \"comparison_samples\": [0.5042158516020236, 0.4904171364148816, 0.49253731343283585]}"
          },
          {
            "name": "ailist/10K/basic/left-adjacent-gap_100%",
            "value": 0.6049469964664311,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00885, 0.0085, 0.00856], \"comparison_samples\": [0.6362329259525521, 0.5952380952380953, 0.6049469964664311]}"
          },
          {
            "name": "ailist/10K/basic/mid-gap1_100%",
            "value": 0.6352517985611511,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00914, 0.00883, 0.00849], \"comparison_samples\": [0.6575539568345324, 0.6316165951359084, 0.6183539694100509]}"
          },
          {
            "name": "ailist/10K/basic/mid-gap2_100%",
            "value": 0.6198406951484432,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00879, 0.00856, 0.00856], \"comparison_samples\": [0.6364952932657494, 0.6109921484653819, 0.6198406951484432]}"
          },
          {
            "name": "ailist/10K/basic/perfect-gap_100%",
            "value": 0.6056826056826057,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00887, 0.00863, 0.00874], \"comparison_samples\": [0.6146916146916146, 0.6051893408134643, 0.5982203969883642]}"
          },
          {
            "name": "ailist/10K/basic/perfect_100%",
            "value": 0.5360934182590233,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01015, 0.00946, 0.0101], \"comparison_samples\": [0.5387473460721868, 0.4952879581151833, 0.5495103373231773]}"
          },
          {
            "name": "ailist/10K/basic/right-adjacent-gap_100%",
            "value": 0.5965156794425087,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00875, 0.00856, 0.0085], \"comparison_samples\": [0.6097560975609757, 0.5940319222761972, 0.6075768406004289]}"
          },
          {
            "name": "ailist/10K/complex/mult_100bin",
            "value": 0.49481621112158347,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00525, 0.00539, 0.00511], \"comparison_samples\": [0.49481621112158347, 0.5089707271010386, 0.4771241830065359]}"
          },
          {
            "name": "awk/10K/basic/3p-partial_100%",
            "value": 77.88530688158711,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.25629, 1.25997, 1.25186], \"comparison_samples\": [76.79034229828851, 80.7155669442665, 77.61066336019839]}"
          },
          {
            "name": "awk/10K/basic/5p-partial_100%",
            "value": 79.39759036144578,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.2521, 1.27035, 1.25122], \"comparison_samples\": [79.39759036144578, 78.75697458152513, 79.44253968253969]}"
          },
          {
            "name": "awk/10K/basic/contained_100%",
            "value": 71.39110604332953,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.24969, 1.25298, 1.2522], \"comparison_samples\": [76.2471018913972, 70.78983050847458, 71.39110604332953]}"
          },
          {
            "name": "awk/10K/basic/enclosed_100%",
            "value": 70.24676784710512,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.24969, 1.25101, 1.24872], \"comparison_samples\": [70.24676784710512, 70.51916572717023, 69.02819237147597]}"
          },
          {
            "name": "awk/10K/basic/left-adjacent-gap_100%",
            "value": 88.97738515901061,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.24835, 1.27181, 1.25903], \"comparison_samples\": [89.74478792235801, 89.06232492997201, 88.97738515901061]}"
          },
          {
            "name": "awk/10K/basic/mid-gap1_100%",
            "value": 90.12230215827338,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.2527, 1.25615, 1.25073], \"comparison_samples\": [90.12230215827338, 89.85336194563664, 91.09468317552803]}"
          },
          {
            "name": "awk/10K/basic/mid-gap2_100%",
            "value": 90.7002172338885,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.25257, 1.25278, 1.25143], \"comparison_samples\": [90.7002172338885, 89.42041399000713, 90.61766835626358]}"
          },
          {
            "name": "awk/10K/basic/perfect-gap_100%",
            "value": 86.72002772002773,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.25137, 1.25908, 1.24781], \"comparison_samples\": [86.72002772002773, 88.2945301542777, 85.40793976728268]}"
          },
          {
            "name": "awk/10K/basic/perfect_100%",
            "value": 66.56794055201699,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.25443, 1.25414, 1.24982], \"comparison_samples\": [66.58333333333333, 65.66178010471205, 67.99891186071817]}"
          },
          {
            "name": "awk/10K/basic/right-adjacent-gap_100%",
            "value": 87.31080139372823,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.25369, 1.25291, 1.25134], \"comparison_samples\": [87.36515679442509, 86.94725884802222, 89.44531808434596]}"
          },
          {
            "name": "awk/10K/complex/mult_100bin",
            "value": 5.908576814326107,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.06269, 0.06547, 0.06207], \"comparison_samples\": [5.908576814326107, 6.182247403210575, 5.795518207282913]}"
          },
          {
            "name": "bedmaps/10K/basic/3p-partial_100%",
            "value": 1.6720396776193431,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02685, 0.02697, 0.02727], \"comparison_samples\": [1.64119804400978, 1.7277386290839205, 1.690638561686299]}"
          },
          {
            "name": "bedmaps/10K/basic/5p-partial_100%",
            "value": 1.6911857958148384,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02667, 0.0266, 0.02807], \"comparison_samples\": [1.6911857958148384, 1.6491010539367639, 1.7822222222222224]}"
          },
          {
            "name": "bedmaps/10K/basic/contained_100%",
            "value": 1.532497149372862,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02688, 0.02663, 0.02712], \"comparison_samples\": [1.6400244051250765, 1.5045197740112994, 1.5461801596351197]}"
          },
          {
            "name": "bedmaps/10K/basic/enclosed_100%",
            "value": 1.5098369870713884,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02756, 0.0263, 0.02686], \"comparison_samples\": [1.5491849353569422, 1.4825253664036078, 1.484798231066888]}"
          },
          {
            "name": "bedmaps/10K/basic/left-adjacent-gap_100%",
            "value": 1.5597173144876326,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02218, 0.02207, 0.02152], \"comparison_samples\": [1.5945363048166785, 1.5455182072829132, 1.5208480565371025]}"
          },
          {
            "name": "bedmaps/10K/basic/mid-gap1_100%",
            "value": 1.566906474820144,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02178, 0.02208, 0.02136], \"comparison_samples\": [1.566906474820144, 1.5793991416309012, 1.555717407137655]}"
          },
          {
            "name": "bedmaps/10K/basic/mid-gap2_100%",
            "value": 1.5778421433743666,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02165, 0.02179, 0.0218], \"comparison_samples\": [1.5677045619116583, 1.5553176302640972, 1.5785662563359886]}"
          },
          {
            "name": "bedmaps/10K/basic/perfect-gap_100%",
            "value": 1.5155925155925156,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02187, 0.02201, 0.02181], \"comparison_samples\": [1.5155925155925156, 1.5434782608695652, 1.4928131416837782]}"
          },
          {
            "name": "bedmaps/10K/basic/perfect_100%",
            "value": 1.4856687898089174,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02789, 0.02799, 0.02856], \"comparison_samples\": [1.4803609341825903, 1.4654450261780105, 1.5538628944504895]}"
          },
          {
            "name": "bedmaps/10K/basic/right-adjacent-gap_100%",
            "value": 1.5303135888501742,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02203, 0.02196, 0.0217], \"comparison_samples\": [1.535191637630662, 1.523941707147814, 1.551107934238742]}"
          },
          {
            "name": "bedmaps/10K/complex/mult_100bin",
            "value": 1.709707822808671,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01814, 0.01893, 0.01786], \"comparison_samples\": [1.709707822808671, 1.7875354107648724, 1.6676003734827265]}"
          },
          {
            "name": "bedops/10K/basic/3p-partial_100%",
            "value": 1.5399876007439555,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02484, 0.02441, 0.025], \"comparison_samples\": [1.5183374083129586, 1.5637411915438821, 1.5499070055796655]}"
          },
          {
            "name": "bedops/10K/basic/5p-partial_100%",
            "value": 1.5668991756499684,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02468, 0.02482, 0.02471], \"comparison_samples\": [1.564996829422955, 1.5387476751394917, 1.5688888888888888]}"
          },
          {
            "name": "bedops/10K/basic/contained_100%",
            "value": 1.4161915621436716,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02481, 0.02484, 0.02518], \"comparison_samples\": [1.5137278828553997, 1.4033898305084747, 1.4355758266818701]}"
          },
          {
            "name": "bedops/10K/basic/enclosed_100%",
            "value": 1.4103428892636314,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02509, 0.02468, 0.02543], \"comparison_samples\": [1.4103428892636314, 1.391206313416009, 1.4057490326147044]}"
          },
          {
            "name": "bedops/10K/basic/left-adjacent-gap_100%",
            "value": 1.5667844522968197,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0222, 0.02217, 0.02201], \"comparison_samples\": [1.5959741193386052, 1.5525210084033614, 1.55547703180212]}"
          },
          {
            "name": "bedops/10K/basic/mid-gap1_100%",
            "value": 1.6035971223021583,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02247, 0.02229, 0.02199], \"comparison_samples\": [1.6165467625899281, 1.5944206008583692, 1.6016023306627822]}"
          },
          {
            "name": "bedops/10K/basic/mid-gap2_100%",
            "value": 1.5865314989138306,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02187, 0.02216, 0.02191], \"comparison_samples\": [1.5836350470673426, 1.5817273376159886, 1.5865314989138306]}"
          },
          {
            "name": "bedops/10K/basic/perfect-gap_100%",
            "value": 1.5176715176715176,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02217, 0.0219, 0.02171], \"comparison_samples\": [1.5363825363825363, 1.5357643758765778, 1.485968514715948]}"
          },
          {
            "name": "bedops/10K/basic/perfect_100%",
            "value": 1.397027600849257,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02643, 0.02632, 0.02631], \"comparison_samples\": [1.4028662420382165, 1.3780104712041885, 1.4314472252448314]}"
          },
          {
            "name": "bedops/10K/basic/right-adjacent-gap_100%",
            "value": 1.575609756097561,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02261, 0.02214, 0.02269], \"comparison_samples\": [1.575609756097561, 1.5364330326162388, 1.6218727662616153]}"
          },
          {
            "name": "bedops/10K/complex/mult_100bin",
            "value": 1.692742695570217,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01796, 0.01808, 0.01764], \"comparison_samples\": [1.692742695570217, 1.7072710103871576, 1.6470588235294117]}"
          },
          {
            "name": "bedtk/10K/basic/3p-partial_100%",
            "value": 0.40855548667079983,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00641, 0.0067, 0.00659], \"comparison_samples\": [0.3918092909535452, 0.42921204356181936, 0.40855548667079983]}"
          },
          {
            "name": "bedtk/10K/basic/5p-partial_100%",
            "value": 0.4134432466708941,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00652, 0.00665, 0.00635], \"comparison_samples\": [0.4134432466708941, 0.41227526348419097, 0.40317460317460313]}"
          },
          {
            "name": "bedtk/10K/basic/contained_100%",
            "value": 0.3694412770809578,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00648, 0.00655, 0.00642], \"comparison_samples\": [0.39536302623550945, 0.3700564971751413, 0.3660205245153934]}"
          },
          {
            "name": "bedtk/10K/basic/enclosed_100%",
            "value": 0.372681281618887,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00663, 0.00696, 0.0063], \"comparison_samples\": [0.372681281618887, 0.39233370913190535, 0.3482587064676617]}"
          },
          {
            "name": "bedtk/10K/basic/left-adjacent-gap_100%",
            "value": 0.39151943462897526,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00554, 0.00559, 0.00546], \"comparison_samples\": [0.39827462257368795, 0.3914565826330533, 0.3858657243816254]}"
          },
          {
            "name": "bedtk/10K/basic/mid-gap1_100%",
            "value": 0.39280575539568346,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00546, 0.00546, 0.00537], \"comparison_samples\": [0.39280575539568346, 0.3905579399141631, 0.3911143481427531]}"
          },
          {
            "name": "bedtk/10K/basic/mid-gap2_100%",
            "value": 0.40260680666183923,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0056, 0.00556, 0.00546], \"comparison_samples\": [0.4055032585083273, 0.396859386152748, 0.3953656770456191]}"
          },
          {
            "name": "bedtk/10K/basic/perfect-gap_100%",
            "value": 0.3825363825363825,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00552, 0.00559, 0.00551], \"comparison_samples\": [0.3825363825363825, 0.39200561009817675, 0.37713894592744696]}"
          },
          {
            "name": "bedtk/10K/basic/perfect_100%",
            "value": 0.41666666666666663,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00785, 0.00865, 0.00747], \"comparison_samples\": [0.41666666666666663, 0.45287958115183247, 0.40642002176278563]}"
          },
          {
            "name": "bedtk/10K/basic/right-adjacent-gap_100%",
            "value": 0.3860627177700348,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00554, 0.00554, 0.00544], \"comparison_samples\": [0.3860627177700348, 0.3844552394170715, 0.3888491779842745]}"
          },
          {
            "name": "bedtk/10K/complex/mult_100bin",
            "value": 0.37417530631479734,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00397, 0.00398, 0.00392], \"comparison_samples\": [0.37417530631479734, 0.37582625118035884, 0.3660130718954248]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/3p-partial_100%",
            "value": 1.078735275883447,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0174, 0.01751, 0.01704], \"comparison_samples\": [1.0635696821515892, 1.1217168481742472, 1.0564166150030998]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/5p-partial_100%",
            "value": 1.105263157894737,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01743, 0.01802, 0.01689], \"comparison_samples\": [1.105263157894737, 1.1171729696218229, 1.0723809523809522]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/contained_100%",
            "value": 0.9561003420752565,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01666, 0.01717, 0.01677], \"comparison_samples\": [1.0164734594264797, 0.9700564971751413, 0.9561003420752565]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/enclosed_100%",
            "value": 0.949971894322653,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01649, 0.0169, 0.0172], \"comparison_samples\": [0.9269252388982575, 0.9526493799323562, 0.9508015478164733]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/left-adjacent-gap_100%",
            "value": 1.1349823321554773,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01563, 0.01606, 0.01698], \"comparison_samples\": [1.1236520488856938, 1.1246498599439778, 1.2]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/mid-gap1_100%",
            "value": 1.1517985611510793,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01601, 0.01591, 0.01613], \"comparison_samples\": [1.1517985611510793, 1.1380543633762519, 1.1747997086671522]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/mid-gap2_100%",
            "value": 1.1426502534395366,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01611, 0.01578, 0.01569], \"comparison_samples\": [1.166545981173063, 1.1263383297644538, 1.1361332367849384]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/perfect-gap_100%",
            "value": 1.1053361053361053,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01585, 0.01595, 0.0166], \"comparison_samples\": [1.0984060984060984, 1.1185133239831697, 1.1362080766598222]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/perfect_100%",
            "value": 0.9670912951167728,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01769, 0.01833, 0.01822], \"comparison_samples\": [0.93895966029724, 0.9596858638743455, 0.9912948857453754]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/right-adjacent-gap_100%",
            "value": 1.1247386759581883,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01614, 0.01624, 0.01599], \"comparison_samples\": [1.1247386759581883, 1.126995142262318, 1.1429592566118656]}"
          },
          {
            "name": "bedtk_sorted/10K/complex/mult_100bin",
            "value": 0.5560791705937794,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0059, 0.00606, 0.00577], \"comparison_samples\": [0.5560791705937794, 0.5722379603399433, 0.5387488328664799]}"
          },
          {
            "name": "bedtools/10K/basic/3p-partial_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01636, 0.01561, 0.01613], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/5p-partial_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01577, 0.01613, 0.01575], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/contained_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01639, 0.0177, 0.01754], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/enclosed_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01779, 0.01774, 0.01809], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/left-adjacent-gap_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01391, 0.01428, 0.01415], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/mid-gap1_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0139, 0.01398, 0.01373], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/mid-gap2_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01381, 0.01401, 0.01381], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/perfect-gap_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01443, 0.01426, 0.01461], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/perfect_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01884, 0.0191, 0.01838], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/right-adjacent-gap_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01435, 0.01441, 0.01399], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/complex/mult_100bin",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01061, 0.01059, 0.01071], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/3p-partial_100%",
            "value": 1.258524488530688,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02012, 0.0203, 0.02081], \"comparison_samples\": [1.2298288508557458, 1.3004484304932733, 1.2901425914445133]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/5p-partial_100%",
            "value": 1.2942295497780598,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02041, 0.02026, 0.02156], \"comparison_samples\": [1.2942295497780598, 1.2560446373217609, 1.3688888888888888]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/contained_100%",
            "value": 1.1892816419612315,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02086, 0.02051, 0.02086], \"comparison_samples\": [1.272727272727273, 1.1587570621468926, 1.1892816419612315]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/enclosed_100%",
            "value": 1.1618887015177066,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02044, 0.02067, 0.02086], \"comparison_samples\": [1.1489600899381676, 1.1651634723788051, 1.1531232725262577]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/left-adjacent-gap_100%",
            "value": 1.3081272084805653,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01883, 0.01833, 0.01851], \"comparison_samples\": [1.3537023723939612, 1.2836134453781514, 1.3081272084805653]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/mid-gap1_100%",
            "value": 1.325179856115108,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01855, 0.01835, 0.01842], \"comparison_samples\": [1.3345323741007196, 1.3125894134477827, 1.3415877640203933]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/mid-gap2_100%",
            "value": 1.3388848660391022,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01839, 0.01849, 0.01933], \"comparison_samples\": [1.3316437364228821, 1.3197715917201998, 1.3997103548153513]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/perfect-gap_100%",
            "value": 1.2889812889812888,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0186, 0.01837, 0.01874], \"comparison_samples\": [1.2889812889812888, 1.288218793828892, 1.2826830937713896]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/perfect_100%",
            "value": 1.1725053078556262,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02209, 0.02237, 0.02203], \"comparison_samples\": [1.1725053078556262, 1.1712041884816755, 1.1985854189336236]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/right-adjacent-gap_100%",
            "value": 1.340766550522648,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01834, 0.01924, 0.0193], \"comparison_samples\": [1.278048780487805, 1.3351839000693964, 1.3795568263045033]}"
          },
          {
            "name": "bedtools_sorted/10K/complex/mult_100bin",
            "value": 1.1281809613572102,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01166, 0.01197, 0.01197], \"comparison_samples\": [1.0989632422243167, 1.130311614730878, 1.1176470588235292]}"
          },
          {
            "name": "genogrove/10K/basic/3p-partial_100%",
            "value": 0.9950402975821451,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01605, 0.01647, 0.01603], \"comparison_samples\": [0.9810513447432762, 1.0550928891736064, 0.9938003719776813]}"
          },
          {
            "name": "genogrove/10K/basic/5p-partial_100%",
            "value": 1.020291693088142,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01609, 0.01642, 0.01581], \"comparison_samples\": [1.020291693088142, 1.0179789212647243, 1.003809523809524]}"
          },
          {
            "name": "genogrove/10K/basic/contained_100%",
            "value": 0.9293044469783351,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01606, 0.01656, 0.0163], \"comparison_samples\": [0.9798657718120807, 0.935593220338983, 0.9293044469783351]}"
          },
          {
            "name": "genogrove/10K/basic/enclosed_100%",
            "value": 0.9151208544125913,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01635, 0.01627, 0.01628], \"comparison_samples\": [0.9190556492411467, 0.9171364148816235, 0.8999447208402432]}"
          },
          {
            "name": "genogrove/10K/basic/left-adjacent-gap_100%",
            "value": 0.9929328621908128,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01405, 0.01526, 0.01398], \"comparison_samples\": [1.0100647016534867, 1.0686274509803921, 0.9879858657243816]}"
          },
          {
            "name": "genogrove/10K/basic/mid-gap1_100%",
            "value": 1.0014388489208634,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01392, 0.01399, 0.01389], \"comparison_samples\": [1.0014388489208634, 1.0007153075822606, 1.0116533139111434]}"
          },
          {
            "name": "genogrove/10K/basic/mid-gap2_100%",
            "value": 1.0130340333091963,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01383, 0.01399, 0.01412], \"comparison_samples\": [1.001448225923244, 0.9985724482512491, 1.0224475018102825]}"
          },
          {
            "name": "genogrove/10K/basic/perfect-gap_100%",
            "value": 0.9812889812889813,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01396, 0.01421, 0.01416], \"comparison_samples\": [0.9674289674289674, 0.9964936886395512, 0.9691991786447639]}"
          },
          {
            "name": "genogrove/10K/basic/perfect_100%",
            "value": 0.9331210191082803,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01732, 0.01758, 0.01805], \"comparison_samples\": [0.9193205944798302, 0.9204188481675393, 0.9820457018498367]}"
          },
          {
            "name": "genogrove/10K/basic/right-adjacent-gap_100%",
            "value": 0.9958188153310105,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01429, 0.01413, 0.01482], \"comparison_samples\": [0.9958188153310105, 0.9805690492713394, 1.0593280914939243]}"
          },
          {
            "name": "genogrove/10K/complex/mult_100bin",
            "value": 1.3628652214891612,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01441, 0.01446, 0.01456], \"comparison_samples\": [1.358152686145146, 1.3654390934844192, 1.3594771241830064]}"
          },
          {
            "name": "gia/10K/basic/3p-partial_100%",
            "value": 0.9237445753254806,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0157, 0.0149, 0.01463], \"comparison_samples\": [0.9596577017114913, 0.9545163356822549, 0.9070055796652202]}"
          },
          {
            "name": "gia/10K/basic/5p-partial_100%",
            "value": 0.9720989220038048,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01537, 0.01481, 0.01533], \"comparison_samples\": [0.9746353836398225, 0.9181649101053938, 0.9733333333333333]}"
          },
          {
            "name": "gia/10K/basic/contained_100%",
            "value": 0.8557582668187002,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01516, 0.01501, 0.01493], \"comparison_samples\": [0.9249542403904821, 0.8480225988700565, 0.8511972633979475]}"
          },
          {
            "name": "gia/10K/basic/enclosed_100%",
            "value": 0.8409218662169757,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01496, 0.01475, 0.01552], \"comparison_samples\": [0.8409218662169757, 0.8314543404735062, 0.8579325594250968]}"
          },
          {
            "name": "gia/10K/basic/left-adjacent-gap_100%",
            "value": 0.9773851590106007,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01317, 0.01383, 0.01391], \"comparison_samples\": [0.9468008626887131, 0.9684873949579833, 0.9830388692579506]}"
          },
          {
            "name": "gia/10K/basic/mid-gap1_100%",
            "value": 0.9582733812949641,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01314, 0.01336, 0.01332], \"comparison_samples\": [0.9453237410071943, 0.955650929899857, 0.970138383102695]}"
          },
          {
            "name": "gia/10K/basic/mid-gap2_100%",
            "value": 0.9833454018826937,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01358, 0.01359, 0.01351], \"comparison_samples\": [0.9833454018826937, 0.9700214132762313, 0.9782766111513396]}"
          },
          {
            "name": "gia/10K/basic/perfect-gap_100%",
            "value": 0.9410949410949411,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01352, 0.01361, 0.01358], \"comparison_samples\": [0.936936936936937, 0.9544179523141655, 0.9295003422313484]}"
          },
          {
            "name": "gia/10K/basic/perfect_100%",
            "value": 0.8625265392781317,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.09658, 0.01591, 0.01625], \"comparison_samples\": [5.126326963906582, 0.8329842931937174, 0.8841131664853101]}"
          },
          {
            "name": "gia/10K/basic/right-adjacent-gap_100%",
            "value": 0.972822299651568,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01344, 0.0142, 0.01396], \"comparison_samples\": [0.9365853658536586, 0.9854267869535046, 0.997855611150822]}"
          },
          {
            "name": "gia/10K/complex/mult_100bin",
            "value": 0.9151743638077285,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00967, 0.00971, 0.00984], \"comparison_samples\": [0.9114043355325165, 0.9169027384324834, 0.9187675070028011]}"
          },
          {
            "name": "giggle/10K/basic/3p-partial_100%",
            "value": 52.33168009919405,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.84411, 0.84686, 0.8425], \"comparison_samples\": [51.595965770171155, 54.251121076233176, 52.23186608803472]}"
          },
          {
            "name": "giggle/10K/basic/5p-partial_100%",
            "value": 53.86239695624604,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.84941, 0.84632, 0.85214], \"comparison_samples\": [53.86239695624604, 52.46869187848729, 54.104126984126985]}"
          },
          {
            "name": "giggle/10K/basic/contained_100%",
            "value": 48.155644241733185,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.8493, 0.84465, 0.84443], \"comparison_samples\": [51.81818181818183, 47.720338983050844, 48.14310148232611]}"
          },
          {
            "name": "giggle/10K/basic/enclosed_100%",
            "value": 47.55424395727937,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.84274, 0.84599, 0.84695], \"comparison_samples\": [47.37155705452501, 47.68827508455468, 46.81868435599779]}"
          },
          {
            "name": "giggle/10K/basic/left-adjacent-gap_100%",
            "value": 3.026148409893993,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04236, 0.04448, 0.04282], \"comparison_samples\": [3.04529115744069, 3.11484593837535, 3.026148409893993]}"
          },
          {
            "name": "giggle/10K/basic/mid-gap1_100%",
            "value": 3.058992805755396,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04252, 0.04796, 0.04242], \"comparison_samples\": [3.058992805755396, 3.430615164520744, 3.0895848506919155]}"
          },
          {
            "name": "giggle/10K/basic/mid-gap2_100%",
            "value": 3.3888486603910213,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04327, 0.04682, 0.0468], \"comparison_samples\": [3.133236784938451, 3.341898643825839, 3.3888486603910213]}"
          },
          {
            "name": "giggle/10K/basic/perfect-gap_100%",
            "value": 3.2266112266112263,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04656, 0.04276, 0.04765], \"comparison_samples\": [3.2266112266112263, 2.9985974754558202, 3.2614647501711156]}"
          },
          {
            "name": "giggle/10K/basic/perfect_100%",
            "value": 45.01751592356688,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.84813, 0.84736, 0.85256], \"comparison_samples\": [45.01751592356688, 44.36439790575916, 46.38520130576713]}"
          },
          {
            "name": "giggle/10K/basic/right-adjacent-gap_100%",
            "value": 3.0766550522648086,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04415, 0.04272, 0.04632], \"comparison_samples\": [3.0766550522648086, 2.9646079111727968, 3.3109363831308074]}"
          },
          {
            "name": "giggle/10K/complex/mult_100bin",
            "value": 173.25824693685203,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.83827, 1.84341, 1.83679], \"comparison_samples\": [173.25824693685203, 174.07082152974505, 171.50233426704014]}"
          },
          {
            "name": "granges/10K/basic/3p-partial_100%",
            "value": 0.9163050216986981,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01492, 0.01467, 0.01478], \"comparison_samples\": [0.9119804400977994, 0.9397821909032672, 0.9163050216986981]}"
          },
          {
            "name": "granges/10K/basic/5p-partial_100%",
            "value": 0.9302473050095118,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01449, 0.01507, 0.01467], \"comparison_samples\": [0.9188332276474318, 0.9342839429634223, 0.9314285714285715]}"
          },
          {
            "name": "granges/10K/basic/contained_100%",
            "value": 0.8363740022805017,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01454, 0.01467, 0.01488], \"comparison_samples\": [0.8871262965222698, 0.8288135593220339, 0.8483466362599772]}"
          },
          {
            "name": "granges/10K/basic/enclosed_100%",
            "value": 0.8302417088251827,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01472, 0.01489, 0.01477], \"comparison_samples\": [0.8274311410905003, 0.8393461104847803, 0.816473189607518]}"
          },
          {
            "name": "granges/10K/basic/left-adjacent-gap_100%",
            "value": 0.6254416961130742,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00885, 0.00927, 0.00866], \"comparison_samples\": [0.6362329259525521, 0.6491596638655462, 0.6120141342756183]}"
          },
          {
            "name": "granges/10K/basic/mid-gap1_100%",
            "value": 0.6280575539568346,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00873, 0.00895, 0.00859], \"comparison_samples\": [0.6280575539568346, 0.6402002861230329, 0.6256372906045158]}"
          },
          {
            "name": "granges/10K/basic/mid-gap2_100%",
            "value": 0.6285300506879073,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00868, 0.00914, 0.00864], \"comparison_samples\": [0.6285300506879073, 0.6523911491791577, 0.6256335988414193]}"
          },
          {
            "name": "granges/10K/basic/perfect-gap_100%",
            "value": 0.6126126126126127,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00882, 0.00895, 0.00884], \"comparison_samples\": [0.6112266112266111, 0.6276297335203366, 0.6050650239561944]}"
          },
          {
            "name": "granges/10K/basic/perfect_100%",
            "value": 0.8710191082802549,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01615, 0.01651, 0.01641], \"comparison_samples\": [0.8572186836518048, 0.8643979057591624, 0.8928182807399347]}"
          },
          {
            "name": "granges/10K/basic/right-adjacent-gap_100%",
            "value": 0.6174216027874564,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00886, 0.00927, 0.0088], \"comparison_samples\": [0.6174216027874564, 0.6433032616238724, 0.6290207290922087]}"
          },
          {
            "name": "granges/10K/complex/mult_100bin",
            "value": 0.9330819981149859,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00965, 0.01036, 0.0099], \"comparison_samples\": [0.9095193213949105, 0.9782813975448535, 0.9243697478991597]}"
          },
          {
            "name": "igd/10K/basic/3p-partial_100%",
            "value": 0.8140111593304402,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0133, 0.01313, 0.01306], \"comparison_samples\": [0.812958435207824, 0.8411274823830877, 0.8096714197148173]}"
          },
          {
            "name": "igd/10K/basic/5p-partial_100%",
            "value": 0.8218135700697526,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01328, 0.01296, 0.01292], \"comparison_samples\": [0.8421052631578948, 0.8034717916924985, 0.8203174603174602]}"
          },
          {
            "name": "igd/10K/basic/contained_100%",
            "value": 0.7599771949828962,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01355, 0.01324, 0.01333], \"comparison_samples\": [0.8267236119585113, 0.7480225988700565, 0.7599771949828962]}"
          },
          {
            "name": "igd/10K/basic/enclosed_100%",
            "value": 0.7830241708825182,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01406, 0.01393, 0.01382], \"comparison_samples\": [0.7903316469926925, 0.7852311161217588, 0.7639579878385849]}"
          },
          {
            "name": "igd/10K/basic/left-adjacent-gap_100%",
            "value": 0.8212014134275619,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01177, 0.01162, 0.01162], \"comparison_samples\": [0.846153846153846, 0.8137254901960784, 0.8212014134275619]}"
          },
          {
            "name": "igd/10K/basic/mid-gap1_100%",
            "value": 0.8424460431654677,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01171, 0.01177, 0.01166], \"comparison_samples\": [0.8424460431654677, 0.8419170243204578, 0.8492352512745813]}"
          },
          {
            "name": "igd/10K/basic/mid-gap2_100%",
            "value": 0.8515568428674873,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01191, 0.01176, 0.01132], \"comparison_samples\": [0.8624185372918176, 0.8394004282655246, 0.8196958725561188]}"
          },
          {
            "name": "igd/10K/basic/perfect-gap_100%",
            "value": 0.8392238392238391,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01211, 0.01194, 0.01344], \"comparison_samples\": [0.8392238392238391, 0.8373071528751752, 0.9199178644763861]}"
          },
          {
            "name": "igd/10K/basic/perfect_100%",
            "value": 0.7914012738853503,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01491, 0.0148, 0.01516], \"comparison_samples\": [0.7914012738853503, 0.774869109947644, 0.82480957562568]}"
          },
          {
            "name": "igd/10K/basic/right-adjacent-gap_100%",
            "value": 0.8118466898954705,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01172, 0.01165, 0.01165], \"comparison_samples\": [0.8167247386759582, 0.808466342817488, 0.8327376697641172]}"
          },
          {
            "name": "igd/10K/complex/mult_100bin",
            "value": 0.7474081055607917,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00803, 0.00775, 0.00793], \"comparison_samples\": [0.756833176248822, 0.7318224740321058, 0.7404295051353874]}"
          },
          {
            "name": "intervaltree/10K/basic/3p-partial_100%",
            "value": 20.046497210167395,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.32104, 0.32948, 0.32335], \"comparison_samples\": [19.623471882640587, 21.106982703395257, 20.046497210167395]}"
          },
          {
            "name": "intervaltree/10K/basic/5p-partial_100%",
            "value": 19.976537729866834,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.31503, 0.32158, 0.31036], \"comparison_samples\": [19.976537729866834, 19.93676379417235, 19.70539682539683]}"
          },
          {
            "name": "intervaltree/10K/basic/contained_100%",
            "value": 15.685290763968071,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.27056, 0.27512, 0.27807], \"comparison_samples\": [16.507626601586335, 15.543502824858756, 15.853477765108323]}"
          },
          {
            "name": "intervaltree/10K/basic/enclosed_100%",
            "value": 20.1186059584036,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.35791, 0.35901, 0.35563], \"comparison_samples\": [20.1186059584036, 20.23731679819617, 19.65892758430072]}"
          },
          {
            "name": "intervaltree/10K/basic/left-adjacent-gap_100%",
            "value": 18.554063604240284,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26254, 0.2629, 0.26071], \"comparison_samples\": [18.874191229331416, 18.410364145658267, 18.424734982332158]}"
          },
          {
            "name": "intervaltree/10K/basic/mid-gap1_100%",
            "value": 18.812230215827338,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.25973, 0.26149, 0.26178], \"comparison_samples\": [18.685611510791368, 18.70457796852647, 19.06627822286963]}"
          },
          {
            "name": "intervaltree/10K/basic/mid-gap2_100%",
            "value": 18.937726285300506,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26131, 0.26163, 0.26153], \"comparison_samples\": [18.92179580014482, 18.674518201284794, 18.937726285300506]}"
          },
          {
            "name": "intervaltree/10K/basic/perfect-gap_100%",
            "value": 18.109494109494108,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.25916, 0.26525, 0.26132], \"comparison_samples\": [17.95980595980596, 18.600981767180926, 17.88637919233402]}"
          },
          {
            "name": "intervaltree/10K/basic/perfect_100%",
            "value": 16.631104033970278,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.31333, 0.31934, 0.31173], \"comparison_samples\": [16.631104033970278, 16.719371727748694, 16.960282916213274]}"
          },
          {
            "name": "intervaltree/10K/basic/right-adjacent-gap_100%",
            "value": 18.18885017421603,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26287, 0.26101, 0.26048], \"comparison_samples\": [18.318466898954703, 18.11311589174185, 18.619013581129376]}"
          },
          {
            "name": "intervaltree/10K/complex/mult_100bin",
            "value": 35.89820923656927,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.38088, 0.38383, 0.37553], \"comparison_samples\": [35.89820923656927, 36.244570349386215, 35.063492063492056]}"
          },
          {
            "name": "polars_bio/10K/basic/3p-partial_100%",
            "value": 16.87042777433354,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.27728, 0.25801, 0.27212], \"comparison_samples\": [16.94865525672372, 16.52850736707239, 16.87042777433354]}"
          },
          {
            "name": "polars_bio/10K/basic/5p-partial_100%",
            "value": 16.666455294863667,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26283, 0.26186, 0.27626], \"comparison_samples\": [16.666455294863667, 16.234345939243646, 17.54031746031746]}"
          },
          {
            "name": "polars_bio/10K/basic/contained_100%",
            "value": 15.091220068415051,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.2647, 0.26161, 0.26924], \"comparison_samples\": [16.150091519219036, 14.780225988700565, 15.350057012542758]}"
          },
          {
            "name": "polars_bio/10K/basic/enclosed_100%",
            "value": 14.76335019673974,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26264, 0.26101, 0.26526], \"comparison_samples\": [14.76335019673974, 14.713077790304398, 14.663349917081261]}"
          },
          {
            "name": "polars_bio/10K/basic/left-adjacent-gap_100%",
            "value": 18.481272084805656,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26151, 0.26005, 0.26181], \"comparison_samples\": [18.800143781452192, 18.21078431372549, 18.502473498233215]}"
          },
          {
            "name": "polars_bio/10K/basic/mid-gap1_100%",
            "value": 19.26978417266187,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26785, 0.25621, 0.27576], \"comparison_samples\": [19.26978417266187, 18.32689556509299, 20.08448652585579]}"
          },
          {
            "name": "polars_bio/10K/basic/mid-gap2_100%",
            "value": 19.224475018102826,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26286, 0.26874, 0.26549], \"comparison_samples\": [19.034033309196236, 19.182012847965737, 19.224475018102826]}"
          },
          {
            "name": "polars_bio/10K/basic/perfect-gap_100%",
            "value": 18.220374220374218,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26363, 0.26292, 0.26059], \"comparison_samples\": [18.26957726957727, 18.43758765778401, 17.836413415468858]}"
          },
          {
            "name": "polars_bio/10K/basic/perfect_100%",
            "value": 14.228237791932061,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.28793, 0.26625, 0.26806], \"comparison_samples\": [15.282908704883228, 13.93979057591623, 14.584330794341676]}"
          },
          {
            "name": "polars_bio/10K/basic/right-adjacent-gap_100%",
            "value": 17.977700348432055,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.25798, 0.25713, 0.26436], \"comparison_samples\": [17.977700348432055, 17.843858431644694, 18.896354538956395]}"
          },
          {
            "name": "polars_bio/10K/complex/mult_100bin",
            "value": 24.75306314797361,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.26976, 0.2572, 0.26263], \"comparison_samples\": [25.42507068803016, 24.287063267233236, 24.5219421101774]}"
          },
          {
            "name": "pygenogrove/10K/basic/3p-partial_100%",
            "value": 10.892746435213889,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.17306, 0.17838, 0.1757], \"comparison_samples\": [10.578239608801956, 11.427290198590647, 10.892746435213889]}"
          },
          {
            "name": "pygenogrove/10K/basic/5p-partial_100%",
            "value": 11.157260621433101,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.17595, 0.17615, 0.17361], \"comparison_samples\": [11.157260621433101, 10.920644761314323, 11.022857142857141]}"
          },
          {
            "name": "pygenogrove/10K/basic/contained_100%",
            "value": 9.97605473204105,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.17323, 0.17681, 0.17498], \"comparison_samples\": [10.569249542403906, 9.989265536723163, 9.97605473204105]}"
          },
          {
            "name": "pygenogrove/10K/basic/enclosed_100%",
            "value": 9.874648679033164,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.17567, 0.17679, 0.17535], \"comparison_samples\": [9.874648679033164, 9.965614430665164, 9.693200663349918]}"
          },
          {
            "name": "pygenogrove/10K/basic/left-adjacent-gap_100%",
            "value": 10.575971731448764,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.14965, 0.14897, 0.15057], \"comparison_samples\": [10.758447160316319, 10.432072829131652, 10.640989399293288]}"
          },
          {
            "name": "pygenogrove/10K/basic/mid-gap1_100%",
            "value": 10.661870503597124,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.15046, 0.14816, 0.1482], \"comparison_samples\": [10.824460431654678, 10.59799713876967, 10.79388201019665]}"
          },
          {
            "name": "pygenogrove/10K/basic/mid-gap2_100%",
            "value": 10.825488776249095,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.14995, 0.1495, 0.14802], \"comparison_samples\": [10.858073859522086, 10.67094932191292, 10.718320057929038]}"
          },
          {
            "name": "pygenogrove/10K/basic/perfect-gap_100%",
            "value": 10.41926541926542,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.15041, 0.14833, 0.15035], \"comparison_samples\": [10.423423423423422, 10.401823281907433, 10.290896646132786]}"
          },
          {
            "name": "pygenogrove/10K/basic/perfect_100%",
            "value": 9.309978768577496,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.2116, 0.17524, 0.1754], \"comparison_samples\": [11.231422505307856, 9.174869109947645, 9.542981501632209]}"
          },
          {
            "name": "pygenogrove/10K/basic/right-adjacent-gap_100%",
            "value": 10.433449477351916,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.14972, 0.15021, 0.14945], \"comparison_samples\": [10.433449477351916, 10.424011103400417, 10.682630450321657]}"
          },
          {
            "name": "pygenogrove/10K/complex/mult_100bin",
            "value": 11.253534401508013,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.11856, 0.12479, 0.1194], \"comparison_samples\": [11.174363807728557, 11.783758262511803, 11.1484593837535]}"
          },
          {
            "name": "tabix/10K/basic/3p-partial_100%",
            "value": 11.419094854308742,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.18419, 0.18418, 0.18515], \"comparison_samples\": [11.258557457212714, 11.798846893017297, 11.478611283323001]}"
          },
          {
            "name": "tabix/10K/basic/5p-partial_100%",
            "value": 11.699429296131896,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.18249, 0.1845, 0.18473], \"comparison_samples\": [11.571972098922005, 11.43831370117793, 11.72888888888889]}"
          },
          {
            "name": "tabix/10K/basic/contained_100%",
            "value": 10.55188141391106,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.18527, 0.18508, 0.18466], \"comparison_samples\": [11.303843807199513, 10.456497175141243, 10.52793614595211]}"
          },
          {
            "name": "tabix/10K/basic/enclosed_100%",
            "value": 10.378302417088252,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.18406, 0.18463, 0.18508], \"comparison_samples\": [10.346261944912872, 10.407553551296505, 10.231066887783307]}"
          },
          {
            "name": "tabix/10K/basic/left-adjacent-gap_100%",
            "value": 12.896819787985867,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.1824, 0.18251, 0.18249], \"comparison_samples\": [13.112868439971244, 12.780812324929974, 12.896819787985867]}"
          },
          {
            "name": "tabix/10K/basic/mid-gap1_100%",
            "value": 13.251079136690647,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.18303, 0.18419, 0.18527], \"comparison_samples\": [13.167625899280576, 13.175250357653791, 13.493809176984705]}"
          },
          {
            "name": "tabix/10K/basic/mid-gap2_100%",
            "value": 13.325850832729907,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.18368, 0.18476, 0.18403], \"comparison_samples\": [13.300506879073136, 13.187723054960744, 13.325850832729907]}"
          },
          {
            "name": "tabix/10K/basic/perfect-gap_100%",
            "value": 12.744975744975743,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.18391, 0.18487, 0.18297], \"comparison_samples\": [12.744975744975743, 12.964235624123422, 12.523613963039015]}"
          },
          {
            "name": "tabix/10K/basic/perfect_100%",
            "value": 9.868365180467093,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.18468, 0.18592, 0.18707], \"comparison_samples\": [9.802547770700638, 9.734031413612566, 10.177910772578889]}"
          },
          {
            "name": "tabix/10K/basic/right-adjacent-gap_100%",
            "value": 12.726829268292683,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.18254, 0.18473, 0.18263], \"comparison_samples\": [12.720557491289199, 12.819569743233867, 13.054324517512507]}"
          },
          {
            "name": "tabix/10K/complex/mult_100bin",
            "value": 1.0141376060320453,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01077, 0.01076, 0.00974], \"comparison_samples\": [1.0150801131008482, 1.016052880075543, 0.9094304388422035]}"
          },
          {
            "name": "ucsc/10K/basic/3p-partial_100%",
            "value": 0.9398636081835091,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01529, 0.0143, 0.01516], \"comparison_samples\": [0.934596577017115, 0.9160794362588084, 0.9398636081835091]}"
          },
          {
            "name": "ucsc/10K/basic/5p-partial_100%",
            "value": 0.9403931515535827,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01483, 0.01454, 0.01494], \"comparison_samples\": [0.9403931515535827, 0.9014259144451334, 0.9485714285714286]}"
          },
          {
            "name": "ucsc/10K/basic/contained_100%",
            "value": 0.8529076396807297,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01497, 0.01496, 0.01473], \"comparison_samples\": [0.9133618059792558, 0.8451977401129943, 0.8397947548460661]}"
          },
          {
            "name": "ucsc/10K/basic/enclosed_100%",
            "value": 0.8369870713884204,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01489, 0.01523, 0.01447], \"comparison_samples\": [0.8369870713884204, 0.858511837655017, 0.7998894416804866]}"
          },
          {
            "name": "ucsc/10K/basic/left-adjacent-gap_100%",
            "value": 0.8890459363957597,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01267, 0.01258, 0.01255], \"comparison_samples\": [0.9108554996405464, 0.8809523809523809, 0.8869257950530036]}"
          },
          {
            "name": "ucsc/10K/basic/mid-gap1_100%",
            "value": 0.9007194244604317,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01252, 0.01268, 0.01228], \"comparison_samples\": [0.9007194244604317, 0.9070100143061517, 0.8943918426802622]}"
          },
          {
            "name": "ucsc/10K/basic/mid-gap2_100%",
            "value": 0.9036929761042722,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01255, 0.01248, 0.01247], \"comparison_samples\": [0.9087617668356265, 0.8907922912205567, 0.9029688631426503]}"
          },
          {
            "name": "ucsc/10K/basic/perfect-gap_100%",
            "value": 0.8856548856548856,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01278, 0.01303, 0.01268], \"comparison_samples\": [0.8856548856548856, 0.9137447405329593, 0.8678986995208762]}"
          },
          {
            "name": "ucsc/10K/basic/perfect_100%",
            "value": 0.8487261146496816,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0166, 0.01507, 0.01599], \"comparison_samples\": [0.881104033970276, 0.7890052356020942, 0.8699673558215452]}"
          },
          {
            "name": "ucsc/10K/basic/right-adjacent-gap_100%",
            "value": 0.8780487804878049,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01265, 0.0126, 0.01247], \"comparison_samples\": [0.8815331010452961, 0.8743927827897294, 0.8913509649749821]}"
          },
          {
            "name": "ucsc/10K/complex/mult_100bin",
            "value": 1.1008482563619226,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01168, 0.01197, 0.01152], \"comparison_samples\": [1.1008482563619226, 1.130311614730878, 1.0756302521008403]}"
          },
          {
            "name": "ailist/all/end-to-end/wall",
            "value": 0.7647040190000212,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.7647040190000212, 0.7598806020000097, 0.7647427129999755], \"comparison_samples\": [0.7647040190000212, 0.7598806020000097, 0.7647427129999755]}"
          },
          {
            "name": "awk/all/end-to-end/wall",
            "value": 13.206823762999989,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [13.206823762999989, 13.274060054000017, 13.205488954000032], \"comparison_samples\": [13.206823762999989, 13.274060054000017, 13.205488954000032]}"
          },
          {
            "name": "bedmaps/all/end-to-end/wall",
            "value": 0.9108437430000436,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.9131260330000259, 0.9108437430000436, 0.9039791719999926], \"comparison_samples\": [0.9131260330000259, 0.9108437430000436, 0.9039791719999926]}"
          },
          {
            "name": "bedops/all/end-to-end/wall",
            "value": 0.8805893910000009,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8819080920000033, 0.8805893910000009, 0.8683346459999939], \"comparison_samples\": [0.8819080920000033, 0.8805893910000009, 0.8683346459999939]}"
          },
          {
            "name": "bedtk/all/end-to-end/wall",
            "value": 0.6888594900000271,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.6888594900000271, 0.7026643910000416, 0.6821654660000149], \"comparison_samples\": [0.6888594900000271, 0.7026643910000416, 0.6821654660000149]}"
          },
          {
            "name": "bedtk_sorted/all/end-to-end/wall",
            "value": 0.8200848710000059,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8200848710000059, 0.825902077999956, 0.8016918330000067], \"comparison_samples\": [0.8200848710000059, 0.825902077999956, 0.8016918330000067]}"
          },
          {
            "name": "bedtools/all/end-to-end/wall",
            "value": 0.7754034309999724,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.7716747309999903, 0.7833757439999545, 0.7754034309999724], \"comparison_samples\": [0.7716747309999903, 0.7833757439999545, 0.7754034309999724]}"
          },
          {
            "name": "bedtools_sorted/all/end-to-end/wall",
            "value": 0.8243295229999603,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8179935370000067, 0.8243295229999603, 0.8346906590000458], \"comparison_samples\": [0.8179935370000067, 0.8243295229999603, 0.8346906590000458]}"
          },
          {
            "name": "genogrove/all/end-to-end/wall",
            "value": 0.7759643680000181,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.7796808900000087, 0.7759643680000181, 0.7586153080000031], \"comparison_samples\": [0.7796808900000087, 0.7759643680000181, 0.7586153080000031]}"
          },
          {
            "name": "gia/all/end-to-end/wall",
            "value": 0.7312983719999693,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.864089454000009, 0.7312983719999693, 0.7234031329999766], \"comparison_samples\": [0.864089454000009, 0.7312983719999693, 0.7234031329999766]}"
          },
          {
            "name": "giggle/all/end-to-end/wall",
            "value": 7.110686921999957,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [7.16618992299999, 7.104575160000024, 7.110686921999957], \"comparison_samples\": [7.16618992299999, 7.104575160000024, 7.110686921999957]}"
          },
          {
            "name": "granges/all/end-to-end/wall",
            "value": 0.7298788699999932,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.927542786999993, 0.7298788699999932, 0.7254903119999767], \"comparison_samples\": [0.927542786999993, 0.7298788699999932, 0.7254903119999767]}"
          },
          {
            "name": "igd/all/end-to-end/wall",
            "value": 0.9735245680000162,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.9796596829999658, 0.9710882550000406, 0.9735245680000162], \"comparison_samples\": [0.9796596829999658, 0.9710882550000406, 0.9735245680000162]}"
          },
          {
            "name": "intervaltree/all/end-to-end/wall",
            "value": 3.8943125479999594,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [3.8788414569999645, 3.9332044540000197, 3.8943125479999594], \"comparison_samples\": [3.8788414569999645, 3.9332044540000197, 3.8943125479999594]}"
          },
          {
            "name": "polars_bio/all/end-to-end/wall",
            "value": 3.5177365819999977,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [3.5177365819999977, 3.4637765260000037, 3.537636770000006], \"comparison_samples\": [3.5177365819999977, 3.4637765260000037, 3.537636770000006]}"
          },
          {
            "name": "pygenogrove/all/end-to-end/wall",
            "value": 2.3352250749999826,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [2.3575620489999665, 2.3352250749999826, 2.3208230220000132], \"comparison_samples\": [2.3575620489999665, 2.3352250749999826, 2.3208230220000132]}"
          },
          {
            "name": "segmeter/all/simulation/wall",
            "value": 1.2942961429999968,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [1.2942961429999968], \"comparison_samples\": [1.2942961429999968]}"
          },
          {
            "name": "tabix/all/end-to-end/wall",
            "value": 2.4972121949999746,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [2.4903845759999967, 2.5026704830000313, 2.4972121949999746], \"comparison_samples\": [2.4903845759999967, 2.5026704830000313, 2.4972121949999746]}"
          },
          {
            "name": "ucsc/all/end-to-end/wall",
            "value": 0.7751552980000156,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.7694653920000292, 0.7785700649999967, 0.7751552980000156], \"comparison_samples\": [0.7694653920000292, 0.7785700649999967, 0.7751552980000156]}"
          }
        ]
      },
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
          "id": "83da71686d8a8bf48255693841f8ceeeba9c386d",
          "message": "Record dirty checkouts and survive empty option arrays in the launcher (#81)\n\n* Record dirty checkouts and survive empty option arrays in the launcher\n\nThe launcher copied the working tree but recorded a clean HEAD as the\nrun's source; uncommitted changes now mark it -dirty, which the\npublication validator rejects. Without optional flags the empty option\narrays were unbound variables under set -u before bash 4.4.\n\n* Document the launcher provenance fixes in the changelog\n\n* Test the dirty mark, count untracked files regardless of git config\n\nThe launcher test now launches from a scratch git repository clean and\nwith an uncommitted file and checks both recorded commits; the status\ncheck forces --untracked-files=normal and treats a failing git status\nas dirty; the validator names the rejected commit.",
          "timestamp": "2026-10-14T23:04:07-05:00",
          "tree_id": "27d4e8503b199319d92c75dd96b4471f3fbdfcba",
          "url": "https://github.com/ylab-hi/segmeter/commit/83da71686d8a8bf48255693841f8ceeeba9c386d"
        },
        "date": 1792037347005,
        "tool": "customSmallerIsBetter",
        "benches": [
          {
            "name": "ailist/10K/basic/3p-partial_100%",
            "value": 0.5061791967044285,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00983, 0.00976, 0.01002], \"comparison_samples\": [0.517913593256059, 0.5025746652935118, 0.5114854517611026]}"
          },
          {
            "name": "ailist/10K/basic/5p-partial_100%",
            "value": 0.49718957588145124,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00973, 0.00971, 0.00994], \"comparison_samples\": [0.5239633817985999, 0.4961676034747062, 0.497248624312156]}"
          },
          {
            "name": "ailist/10K/basic/contained_100%",
            "value": 0.47596382674916704,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01, 0.00989, 0.01023], \"comparison_samples\": [0.48899755501222497, 0.4641013608634444, 0.48691099476439786]}"
          },
          {
            "name": "ailist/10K/basic/enclosed_100%",
            "value": 0.4479213232007152,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01002, 0.00975, 0.01003], \"comparison_samples\": [0.4726415094339622, 0.43585158694680376, 0.43326133909287257]}"
          },
          {
            "name": "ailist/10K/basic/left-adjacent-gap_100%",
            "value": 0.5572222222222223,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01005, 0.00988, 0.01003], \"comparison_samples\": [0.5894428152492669, 0.5485841199333703, 0.5572222222222223]}"
          },
          {
            "name": "ailist/10K/basic/mid-gap1_100%",
            "value": 0.576524741081703,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01002, 0.00989, 0.01003], \"comparison_samples\": [0.610231425091352, 0.5690448791714614, 0.5744558991981673]}"
          },
          {
            "name": "ailist/10K/basic/mid-gap2_100%",
            "value": 0.5592739648326717,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00984, 0.00986, 0.00988], \"comparison_samples\": [0.5902819436112777, 0.5592739648326717, 0.5449531163816878]}"
          },
          {
            "name": "ailist/10K/basic/perfect-gap_100%",
            "value": 0.5299465240641711,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00991, 0.00982, 0.01001], \"comparison_samples\": [0.5675830469644902, 0.5215082315454063, 0.5352941176470588]}"
          },
          {
            "name": "ailist/10K/basic/perfect_100%",
            "value": 0.4645852749301025,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01009, 0.00965, 0.00997], \"comparison_samples\": [0.4823135755258126, 0.4404381560931082, 0.4645852749301025]}"
          },
          {
            "name": "ailist/10K/basic/right-adjacent-gap_100%",
            "value": 0.550764192139738,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01015, 0.01009, 0.00998], \"comparison_samples\": [0.6074207061639737, 0.550764192139738, 0.5294429708222812]}"
          },
          {
            "name": "ailist/10K/complex/mult_100bin",
            "value": 0.4930340557275542,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00637, 0.00626, 0.00653], \"comparison_samples\": [0.4930340557275542, 0.48982785602503914, 0.44603825136612024]}"
          },
          {
            "name": "awk/10K/basic/3p-partial_100%",
            "value": 76.24407826982493,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.47507, 1.53286, 1.48066], \"comparison_samples\": [77.71707060063225, 78.93202883625128, 75.58244002041859]}"
          },
          {
            "name": "awk/10K/basic/5p-partial_100%",
            "value": 76.03423607562596,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.48642, 1.48799, 1.53592], \"comparison_samples\": [80.04415724286484, 76.03423607562596, 76.8344172086043]}"
          },
          {
            "name": "awk/10K/basic/contained_100%",
            "value": 71.70347453593527,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.49746, 1.50649, 1.5233], \"comparison_samples\": [73.22542787286064, 70.69404035664009, 72.50356972870063]}"
          },
          {
            "name": "awk/10K/basic/enclosed_100%",
            "value": 67.59588734912829,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.53152, 1.51212, 1.49927], \"comparison_samples\": [72.24150943396226, 67.59588734912829, 64.76328293736502]}"
          },
          {
            "name": "awk/10K/basic/left-adjacent-gap_100%",
            "value": 83.62666666666667,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.50106, 1.54, 1.50528], \"comparison_samples\": [88.03870967741936, 85.50805108273181, 83.62666666666667]}"
          },
          {
            "name": "awk/10K/basic/mid-gap1_100%",
            "value": 85.95914844649022,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.48977, 1.50983, 1.49397], \"comparison_samples\": [90.72898903775884, 86.8716915995397, 85.56529209621993]}"
          },
          {
            "name": "awk/10K/basic/mid-gap2_100%",
            "value": 85.33749290981282,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.49992, 1.5076, 1.5045], \"comparison_samples\": [89.97720455908818, 85.51332955190017, 82.98400441257584]}"
          },
          {
            "name": "awk/10K/basic/perfect-gap_100%",
            "value": 80.6,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.49174, 1.53168, 1.50722], \"comparison_samples\": [85.43757159221077, 81.3425385023898, 80.6]}"
          },
          {
            "name": "awk/10K/basic/perfect_100%",
            "value": 69.72833178005591,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.49637, 1.51623, 1.48873], \"comparison_samples\": [71.52820267686424, 69.20264719306253, 69.37232059645854]}"
          },
          {
            "name": "awk/10K/basic/right-adjacent-gap_100%",
            "value": 81.94978165938865,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.50132, 1.49604, 1.50195], \"comparison_samples\": [89.84560143626571, 81.66157205240175, 79.6790450928382]}"
          },
          {
            "name": "awk/10K/complex/mult_100bin",
            "value": 5.880030959752322,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.07548, 0.0772, 0.07597], \"comparison_samples\": [5.842105263157896, 6.040688575899844, 5.189207650273223]}"
          },
          {
            "name": "bedmaps/10K/basic/3p-partial_100%",
            "value": 1.4351184346035015,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02759, 0.02787, 0.02856], \"comparison_samples\": [1.4536354056902, 1.4351184346035015, 1.457886676875957]}"
          },
          {
            "name": "bedmaps/10K/basic/5p-partial_100%",
            "value": 1.4154317833418497,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02765, 0.0277, 0.02801], \"comparison_samples\": [1.4889606892837912, 1.4154317833418497, 1.40120060030015]}"
          },
          {
            "name": "bedmaps/10K/basic/contained_100%",
            "value": 1.3574488338886244,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02773, 0.02852, 0.02865], \"comparison_samples\": [1.3559902200488998, 1.338338808071328, 1.3636363636363635]}"
          },
          {
            "name": "bedmaps/10K/basic/enclosed_100%",
            "value": 1.2717925793473401,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.028, 0.0285, 0.02845], \"comparison_samples\": [1.320754716981132, 1.2740277156906572, 1.2289416846652268]}"
          },
          {
            "name": "bedmaps/10K/basic/left-adjacent-gap_100%",
            "value": 1.2855555555555558,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02302, 0.02314, 0.02363], \"comparison_samples\": [1.3501466275659824, 1.2848417545807884, 1.312777777777778]}"
          },
          {
            "name": "bedmaps/10K/basic/mid-gap1_100%",
            "value": 1.3176064441887227,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02283, 0.0229, 0.02325], \"comparison_samples\": [1.3903775883069427, 1.3176064441887227, 1.331615120274914]}"
          },
          {
            "name": "bedmaps/10K/basic/mid-gap2_100%",
            "value": 1.301191151446398,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02257, 0.02294, 0.02312], \"comparison_samples\": [1.3539292141571684, 1.301191151446398, 1.2752344180915611]}"
          },
          {
            "name": "bedmaps/10K/basic/perfect-gap_100%",
            "value": 1.2550802139037434,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02347, 0.02298, 0.02582], \"comparison_samples\": [1.3442153493699887, 1.2203929899097186, 1.3807486631016042]}"
          },
          {
            "name": "bedmaps/10K/basic/perfect_100%",
            "value": 1.2954333643988816,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0278, 0.02758, 0.0282], \"comparison_samples\": [1.32887189292543, 1.258785942492013, 1.3140726933830382]}"
          },
          {
            "name": "bedmaps/10K/basic/right-adjacent-gap_100%",
            "value": 1.2936681222707422,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0232, 0.0237, 0.02411], \"comparison_samples\": [1.3883901855176541, 1.2936681222707422, 1.2790450928381965]}"
          },
          {
            "name": "bedmaps/10K/complex/mult_100bin",
            "value": 1.4698142414860682,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01818, 0.01904, 0.01899], \"comparison_samples\": [1.407120743034056, 1.4898278560250393, 1.2971311475409837]}"
          },
          {
            "name": "bedops/10K/basic/3p-partial_100%",
            "value": 1.3305870236869206,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02597, 0.02559, 0.02584], \"comparison_samples\": [1.368282402528978, 1.3177136972193615, 1.3190403266972945]}"
          },
          {
            "name": "bedops/10K/basic/5p-partial_100%",
            "value": 1.312212570260603,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02572, 0.02532, 0.02568], \"comparison_samples\": [1.3850296176628971, 1.2938170669391926, 1.2846423211605804]}"
          },
          {
            "name": "bedops/10K/basic/contained_100%",
            "value": 1.2498810090433126,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02655, 0.02573, 0.02626], \"comparison_samples\": [1.2982885085574574, 1.2074143594556546, 1.2498810090433126]}"
          },
          {
            "name": "bedops/10K/basic/enclosed_100%",
            "value": 1.174787662047385,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0263, 0.02578, 0.02628], \"comparison_samples\": [1.2405660377358492, 1.1524362986142154, 1.1352051835853132]}"
          },
          {
            "name": "bedops/10K/basic/left-adjacent-gap_100%",
            "value": 1.2855555555555558,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02354, 0.02313, 0.02314], \"comparison_samples\": [1.3806451612903226, 1.2842865074958356, 1.2855555555555558]}"
          },
          {
            "name": "bedops/10K/basic/mid-gap1_100%",
            "value": 1.3452243958573074,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02313, 0.02354, 0.02338], \"comparison_samples\": [1.4086479902557856, 1.3544303797468353, 1.339060710194731]}"
          },
          {
            "name": "bedops/10K/basic/mid-gap2_100%",
            "value": 1.314804310833806,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02291, 0.02475, 0.02318], \"comparison_samples\": [1.3743251349730052, 1.4038570618264323, 1.2785438499724213]}"
          },
          {
            "name": "bedops/10K/basic/perfect-gap_100%",
            "value": 1.23475935828877,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02332, 0.02309, 0.02279], \"comparison_samples\": [1.3356242840778925, 1.22623473181094, 1.218716577540107]}"
          },
          {
            "name": "bedops/10K/basic/perfect_100%",
            "value": 1.2073625349487418,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02628, 0.02554, 0.02591], \"comparison_samples\": [1.2562141491395793, 1.1656777727065268, 1.2073625349487418]}"
          },
          {
            "name": "bedops/10K/basic/right-adjacent-gap_100%",
            "value": 1.2778384279475983,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02356, 0.02341, 0.02326], \"comparison_samples\": [1.409934171154997, 1.2778384279475983, 1.23395225464191]}"
          },
          {
            "name": "bedops/10K/complex/mult_100bin",
            "value": 1.455108359133127,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01874, 0.01895, 0.0188], \"comparison_samples\": [1.4504643962848298, 1.4827856025039126, 1.284153005464481]}"
          },
          {
            "name": "bedtk/10K/basic/3p-partial_100%",
            "value": 0.3805355303810505,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00739, 0.00739, 0.00737], \"comparison_samples\": [0.3893572181243414, 0.3805355303810505, 0.3762123532414497]}"
          },
          {
            "name": "bedtk/10K/basic/5p-partial_100%",
            "value": 0.3760858456821666,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00755, 0.00736, 0.00732], \"comparison_samples\": [0.40656973613354874, 0.3760858456821666, 0.3661830915457729]}"
          },
          {
            "name": "bedtk/10K/basic/contained_100%",
            "value": 0.3541170871013803,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00737, 0.00765, 0.00744], \"comparison_samples\": [0.36039119804400976, 0.3589863913655561, 0.3541170871013803]}"
          },
          {
            "name": "bedtk/10K/basic/enclosed_100%",
            "value": 0.33705856057219485,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0076, 0.00754, 0.00734], \"comparison_samples\": [0.3584905660377358, 0.33705856057219485, 0.31706263498920084]}"
          },
          {
            "name": "bedtk/10K/basic/left-adjacent-gap_100%",
            "value": 0.3466666666666667,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00618, 0.00635, 0.00624], \"comparison_samples\": [0.3624633431085044, 0.3525818989450305, 0.3466666666666667]}"
          },
          {
            "name": "bedtk/10K/basic/mid-gap1_100%",
            "value": 0.356156501726122,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00601, 0.00619, 0.00633], \"comparison_samples\": [0.3660170523751522, 0.356156501726122, 0.36254295532646047]}"
          },
          {
            "name": "bedtk/10K/basic/mid-gap2_100%",
            "value": 0.3499716392512762,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00611, 0.00646, 0.00617], \"comparison_samples\": [0.36652669466106774, 0.3664208735110607, 0.3403199117484832]}"
          },
          {
            "name": "bedtk/10K/basic/perfect-gap_100%",
            "value": 0.33155080213903737,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00647, 0.00616, 0.0062], \"comparison_samples\": [0.37056128293241697, 0.3271375464684015, 0.33155080213903737]}"
          },
          {
            "name": "bedtk/10K/basic/perfect_100%",
            "value": 0.3406337371854613,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00731, 0.00748, 0.00729], \"comparison_samples\": [0.34942638623326955, 0.34139662254678227, 0.33970177073625346]}"
          },
          {
            "name": "bedtk/10K/basic/right-adjacent-gap_100%",
            "value": 0.33406113537117904,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00609, 0.00633, 0.00612], \"comparison_samples\": [0.3644524236983842, 0.3455240174672489, 0.3246684350132626]}"
          },
          {
            "name": "bedtk/10K/complex/mult_100bin",
            "value": 0.38544891640866874,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00511, 0.00498, 0.00493], \"comparison_samples\": [0.3955108359133127, 0.38967136150234744, 0.33674863387978143]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/3p-partial_100%",
            "value": 0.964984552008239,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01899, 0.01874, 0.01842], \"comparison_samples\": [1.000526870389884, 0.964984552008239, 0.9402756508422664]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/5p-partial_100%",
            "value": 0.9494123658661215,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01858, 0.01867, 0.01853], \"comparison_samples\": [1.0005385029617662, 0.9540112416964741, 0.9269634817408705]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/contained_100%",
            "value": 0.9067110899571632,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01905, 0.01921, 0.01876], \"comparison_samples\": [0.9315403422982885, 0.9014547160957298, 0.8929081389814373]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/enclosed_100%",
            "value": 0.8502458649977648,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01902, 0.02044, 0.01875], \"comparison_samples\": [0.8971698113207547, 0.913723737147966, 0.8099352051835853]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/left-adjacent-gap_100%",
            "value": 1.005,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01769, 0.01809, 0.0185], \"comparison_samples\": [1.0375366568914957, 1.0044419766796222, 1.027777777777778]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/mid-gap1_100%",
            "value": 1.0287686996547756,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01751, 0.01825, 0.01788], \"comparison_samples\": [1.0663824604141292, 1.0500575373993095, 1.0240549828178693]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/mid-gap2_100%",
            "value": 1.0266591038003405,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01723, 0.0181, 0.01854], \"comparison_samples\": [1.0335932813437312, 1.0266591038003405, 1.02261445118588]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/perfect-gap_100%",
            "value": 0.9604278074866309,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01796, 0.01793, 0.01802], \"comparison_samples\": [1.0286368843069875, 0.9522039298990973, 0.9636363636363636]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/perfect_100%",
            "value": 0.8723205964585276,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01885, 0.01872, 0.01865], \"comparison_samples\": [0.9010516252390056, 0.8544043815609311, 0.8690587138863001]}"
          },
          {
            "name": "bedtk_sorted/10K/basic/right-adjacent-gap_100%",
            "value": 1.0049126637554584,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01774, 0.01841, 0.01859], \"comparison_samples\": [1.06163973668462, 1.0049126637554584, 0.9862068965517241]}"
          },
          {
            "name": "bedtk_sorted/10K/complex/mult_100bin",
            "value": 0.6083591331269351,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.00755, 0.00786, 0.00816], \"comparison_samples\": [0.5843653250773995, 0.6150234741784039, 0.5573770491803279]}"
          },
          {
            "name": "bedtools/10K/basic/3p-partial_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01898, 0.01942, 0.01959], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/5p-partial_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01857, 0.01957, 0.01999], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/contained_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02045, 0.02131, 0.02101], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/enclosed_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0212, 0.02237, 0.02315], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/left-adjacent-gap_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01705, 0.01801, 0.018], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/mid-gap1_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01642, 0.01738, 0.01746], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/mid-gap2_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01667, 0.01763, 0.01813], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/perfect-gap_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01746, 0.01883, 0.0187], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/perfect_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02092, 0.02191, 0.02146], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/basic/right-adjacent-gap_100%",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01671, 0.01832, 0.01885], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools/10K/complex/mult_100bin",
            "value": 1,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01292, 0.01278, 0.01464], \"comparison_samples\": [1.0, 1.0, 1.0]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/3p-partial_100%",
            "value": 1.2173017507723998,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02284, 0.02364, 0.02494], \"comparison_samples\": [1.2033719704952581, 1.2173017507723998, 1.2730985196528841]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/5p-partial_100%",
            "value": 1.2100153295861011,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02296, 0.02368, 0.02393], \"comparison_samples\": [1.2364028002154013, 1.2100153295861011, 1.1970985492746373]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/contained_100%",
            "value": 1.1465968586387434,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02348, 0.02409, 0.02457], \"comparison_samples\": [1.1481662591687043, 1.1304551853589864, 1.1694431223227035]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/enclosed_100%",
            "value": 1.0433616450603487,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02321, 0.02397, 0.02334], \"comparison_samples\": [1.094811320754717, 1.0715243629861422, 1.008207343412527]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/left-adjacent-gap_100%",
            "value": 1.1866666666666668,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02098, 0.02136, 0.02186], \"comparison_samples\": [1.2304985337243401, 1.1860077734591892, 1.2144444444444447]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/mid-gap1_100%",
            "value": 1.2410817031070196,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02088, 0.02158, 0.02157], \"comparison_samples\": [1.271619975639464, 1.241657077100115, 1.2353951890034365]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/mid-gap2_100%",
            "value": 1.2144072603516731,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02141, 0.02163, 0.02122], \"comparison_samples\": [1.284343131373725, 1.2268859897901305, 1.1704357418643132]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/perfect-gap_100%",
            "value": 1.1310160427807485,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02094, 0.02115, 0.02154], \"comparison_samples\": [1.1993127147766323, 1.123207647371216, 1.1518716577540107]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/perfect_100%",
            "value": 1.0899347623485556,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02339, 0.0236, 0.02327], \"comparison_samples\": [1.1180688336520077, 1.0771337288909175, 1.0843429636533084]}"
          },
          {
            "name": "bedtools_sorted/10K/basic/right-adjacent-gap_100%",
            "value": 1.1528384279475983,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02093, 0.02112, 0.02136], \"comparison_samples\": [1.2525433871932976, 1.1528384279475983, 1.1331564986737401]}"
          },
          {
            "name": "bedtools_sorted/10K/complex/mult_100bin",
            "value": 1.1857585139318887,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01532, 0.01519, 0.01542], \"comparison_samples\": [1.1857585139318887, 1.1885758998435056, 1.0532786885245902]}"
          },
          {
            "name": "genogrove/10K/basic/3p-partial_100%",
            "value": 0.9984552008238929,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01893, 0.01939, 0.02117], \"comparison_samples\": [0.9973656480505795, 0.9984552008238929, 1.0806533945890762]}"
          },
          {
            "name": "genogrove/10K/basic/5p-partial_100%",
            "value": 0.9473684210526316,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01854, 0.01844, 0.01947], \"comparison_samples\": [0.9983844911147012, 0.9422585590189065, 0.9739869934967484]}"
          },
          {
            "name": "genogrove/10K/basic/contained_100%",
            "value": 0.9157544026653974,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01924, 0.01902, 0.01926], \"comparison_samples\": [0.9408312958435208, 0.8925387142186767, 0.9167063303188957]}"
          },
          {
            "name": "genogrove/10K/basic/enclosed_100%",
            "value": 0.8681269557443003,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01906, 0.01942, 0.02001], \"comparison_samples\": [0.8990566037735849, 0.8681269557443003, 0.8643628509719222]}"
          },
          {
            "name": "genogrove/10K/basic/left-adjacent-gap_100%",
            "value": 0.9127777777777779,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01643, 0.01654, 0.01643], \"comparison_samples\": [0.9636363636363637, 0.9183786785119377, 0.9127777777777779]}"
          },
          {
            "name": "genogrove/10K/basic/mid-gap1_100%",
            "value": 0.9533947065592636,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01657, 0.01624, 0.02082], \"comparison_samples\": [1.0091352009744214, 0.9344073647871117, 1.1924398625429555]}"
          },
          {
            "name": "genogrove/10K/basic/mid-gap2_100%",
            "value": 0.9285309132161088,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01661, 0.01631, 0.01637], \"comparison_samples\": [0.9964007198560287, 0.925127623369257, 0.90292333149476]}"
          },
          {
            "name": "genogrove/10K/basic/perfect-gap_100%",
            "value": 0.8684491978609625,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01624, 0.0162, 0.01643], \"comparison_samples\": [0.9301260022909508, 0.8603292618162507, 0.8786096256684491]}"
          },
          {
            "name": "genogrove/10K/basic/perfect_100%",
            "value": 0.8783783783783783,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01885, 0.01848, 0.01982], \"comparison_samples\": [0.9010516252390056, 0.8434504792332269, 0.9235787511649581]}"
          },
          {
            "name": "genogrove/10K/basic/right-adjacent-gap_100%",
            "value": 0.9066593886462883,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01661, 0.0167, 0.01651], \"comparison_samples\": [0.9940155595451826, 0.9115720524017468, 0.8758620689655173]}"
          },
          {
            "name": "genogrove/10K/complex/mult_100bin",
            "value": 1.3382352941176472,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01776, 0.01711, 0.01729], \"comparison_samples\": [1.3746130030959753, 1.3388106416275432, 1.1810109289617485]}"
          },
          {
            "name": "gia/10K/basic/3p-partial_100%",
            "value": 0.9248197734294542,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01782, 0.01796, 0.0192], \"comparison_samples\": [0.9388830347734457, 0.9248197734294542, 0.9800918836140887]}"
          },
          {
            "name": "gia/10K/basic/5p-partial_100%",
            "value": 0.9611650485436893,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01827, 0.01881, 0.01932], \"comparison_samples\": [0.9838449111470114, 0.9611650485436893, 0.9664832416208103]}"
          },
          {
            "name": "gia/10K/basic/contained_100%",
            "value": 0.8491194669205141,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01738, 0.01981, 0.01784], \"comparison_samples\": [0.8498777506112469, 0.9296105114969498, 0.8491194669205141]}"
          },
          {
            "name": "gia/10K/basic/enclosed_100%",
            "value": 0.8985248100134108,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01923, 0.02014, 0.0201], \"comparison_samples\": [0.9070754716981132, 0.9003129190880644, 0.8682505399568035]}"
          },
          {
            "name": "gia/10K/basic/left-adjacent-gap_100%",
            "value": 0.9527777777777777,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01715, 0.01792, 0.01564], \"comparison_samples\": [1.0058651026392962, 0.9950027762354245, 0.868888888888889]}"
          },
          {
            "name": "gia/10K/basic/mid-gap1_100%",
            "value": 0.9654775604142692,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01678, 0.01809, 0.01573], \"comparison_samples\": [1.0219244823386113, 1.0408515535097813, 0.9009163802978236]}"
          },
          {
            "name": "gia/10K/basic/mid-gap2_100%",
            "value": 1.0005672149744753,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0182, 0.01764, 0.01522], \"comparison_samples\": [1.0917816436712657, 1.0005672149744753, 0.839492553778268]}"
          },
          {
            "name": "gia/10K/basic/perfect-gap_100%",
            "value": 0.843850267379679,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01596, 0.0154, 0.01578], \"comparison_samples\": [0.9140893470790378, 0.8178438661710038, 0.843850267379679]}"
          },
          {
            "name": "gia/10K/basic/perfect_100%",
            "value": 0.9147250698974838,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0357, 0.01963, 0.01908], \"comparison_samples\": [1.7065009560229445, 0.8959379278868098, 0.8890959925442684]}"
          },
          {
            "name": "gia/10K/basic/right-adjacent-gap_100%",
            "value": 0.8995633187772927,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01648, 0.01769, 0.01571], \"comparison_samples\": [0.98623578695392, 0.9656113537117905, 0.83342175066313]}"
          },
          {
            "name": "gia/10K/complex/mult_100bin",
            "value": 0.9891640866873065,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01278, 0.0127, 0.01299], \"comparison_samples\": [0.9891640866873065, 0.9937402190923318, 0.8872950819672131]}"
          },
          {
            "name": "giggle/10K/basic/3p-partial_100%",
            "value": 53.07003089598353,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.03062, 1.03489, 1.0293], \"comparison_samples\": [54.30031612223394, 53.289907312049436, 52.54211332312405]}"
          },
          {
            "name": "giggle/10K/basic/5p-partial_100%",
            "value": 52.35922330097087,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.02467, 1.01858, 1.02539], \"comparison_samples\": [55.178782983306405, 52.048032703117016, 51.29514757378689]}"
          },
          {
            "name": "giggle/10K/basic/contained_100%",
            "value": 48.43788672060923,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.01768, 1.03105, 1.01577], \"comparison_samples\": [49.764303178484106, 48.38338808071328, 48.34697762970014]}"
          },
          {
            "name": "giggle/10K/basic/enclosed_100%",
            "value": 45.97228430934287,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.02059, 1.05715, 1.0284], \"comparison_samples\": [48.14103773584906, 47.25748770675011, 44.423326133909285]}"
          },
          {
            "name": "giggle/10K/basic/left-adjacent-gap_100%",
            "value": 2.835555555555556,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.05027, 0.05189, 0.05104], \"comparison_samples\": [2.948387096774194, 2.8811771238200996, 2.835555555555556]}"
          },
          {
            "name": "giggle/10K/basic/mid-gap1_100%",
            "value": 2.9850402761795167,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.05253, 0.05117, 0.05188], \"comparison_samples\": [3.1991473812423874, 2.944188722669735, 2.9713631156930127]}"
          },
          {
            "name": "giggle/10K/basic/mid-gap2_100%",
            "value": 2.920022688598979,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.04972, 0.0524, 0.05148], \"comparison_samples\": [2.982603479304139, 2.972206466250709, 2.839492553778268]}"
          },
          {
            "name": "giggle/10K/basic/perfect-gap_100%",
            "value": 2.781283422459893,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.05286, 0.04919, 0.05201], \"comparison_samples\": [3.027491408934708, 2.6123207647371216, 2.781283422459893]}"
          },
          {
            "name": "giggle/10K/basic/perfect_100%",
            "value": 47.554520037278664,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [1.02085, 1.01863, 1.02052], \"comparison_samples\": [48.79780114722753, 46.49155636695573, 47.554520037278664]}"
          },
          {
            "name": "giggle/10K/basic/right-adjacent-gap_100%",
            "value": 2.767467248908297,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.05887, 0.05062, 0.0507], \"comparison_samples\": [3.5230400957510475, 2.763100436681223, 2.6896551724137936]}"
          },
          {
            "name": "giggle/10K/complex/mult_100bin",
            "value": 173.22987616099073,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [2.23833, 2.23628, 2.23813], \"comparison_samples\": [173.2453560371517, 174.9827856025039, 152.87773224043715]}"
          },
          {
            "name": "granges/10K/basic/3p-partial_100%",
            "value": 0.9186405767250259,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01799, 0.01784, 0.01757], \"comparison_samples\": [0.9478398314014752, 0.9186405767250259, 0.8968861664114344]}"
          },
          {
            "name": "granges/10K/basic/5p-partial_100%",
            "value": 0.9131323454266734,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01741, 0.01787, 0.01816], \"comparison_samples\": [0.9375336564351103, 0.9131323454266734, 0.9084542271135567]}"
          },
          {
            "name": "granges/10K/basic/contained_100%",
            "value": 0.8519752498810089,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01788, 0.01801, 0.0179], \"comparison_samples\": [0.8743276283618582, 0.8451431252932896, 0.8519752498810089]}"
          },
          {
            "name": "granges/10K/basic/enclosed_100%",
            "value": 0.8211890925346447,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01837, 0.01847, 0.01747], \"comparison_samples\": [0.8665094339622642, 0.8256593652212785, 0.7546436285097192]}"
          },
          {
            "name": "granges/10K/basic/left-adjacent-gap_100%",
            "value": 0.5972222222222222,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01092, 0.01075, 0.01062], \"comparison_samples\": [0.6404692082111437, 0.5968906163242642, 0.59]}"
          },
          {
            "name": "granges/10K/basic/mid-gap1_100%",
            "value": 0.6219792865362486,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01094, 0.01081, 0.01055], \"comparison_samples\": [0.6662606577344702, 0.6219792865362486, 0.6042382588774342]}"
          },
          {
            "name": "granges/10K/basic/mid-gap2_100%",
            "value": 0.6108905275099262,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01114, 0.01077, 0.01048], \"comparison_samples\": [0.6682663467306539, 0.6108905275099262, 0.5780474351902923]}"
          },
          {
            "name": "granges/10K/basic/perfect-gap_100%",
            "value": 0.5652406417112299,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01049, 0.01099, 0.01057], \"comparison_samples\": [0.6008018327605956, 0.5836431226765799, 0.5652406417112299]}"
          },
          {
            "name": "granges/10K/basic/perfect_100%",
            "value": 0.8191985088536812,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01758, 0.01839, 0.01745], \"comparison_samples\": [0.8403441682600381, 0.8393427658603377, 0.8131407269338304]}"
          },
          {
            "name": "granges/10K/basic/right-adjacent-gap_100%",
            "value": 0.5846069868995634,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01071, 0.01087, 0.01054], \"comparison_samples\": [0.6409335727109516, 0.5933406113537117, 0.5591511936339524]}"
          },
          {
            "name": "granges/10K/complex/mult_100bin",
            "value": 1.0278637770897834,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01286, 0.01328, 0.01358], \"comparison_samples\": [0.9953560371517028, 1.0391236306729266, 0.9275956284153005]}"
          },
          {
            "name": "igd/10K/basic/3p-partial_100%",
            "value": 1.115859938208033,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02155, 0.02181, 0.02167], \"comparison_samples\": [1.1354056902002108, 1.123069001029866, 1.1061766207248596]}"
          },
          {
            "name": "igd/10K/basic/5p-partial_100%",
            "value": 1.0991313234542668,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02151, 0.02165, 0.0214], \"comparison_samples\": [1.1583198707592892, 1.1062851303014818, 1.070535267633817]}"
          },
          {
            "name": "igd/10K/basic/contained_100%",
            "value": 1.0471204188481675,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.022, 0.02204, 0.02177], \"comparison_samples\": [1.075794621026895, 1.0342562177381511, 1.0361732508329367]}"
          },
          {
            "name": "igd/10K/basic/enclosed_100%",
            "value": 1.0335270451497542,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02287, 0.0232, 0.02312], \"comparison_samples\": [1.0787735849056603, 1.037103263299061, 0.9987041036717064]}"
          },
          {
            "name": "igd/10K/basic/left-adjacent-gap_100%",
            "value": 1.1111111111111112,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02052, 0.02, 0.01983], \"comparison_samples\": [1.2035190615835778, 1.110494169905608, 1.1016666666666668]}"
          },
          {
            "name": "igd/10K/basic/mid-gap1_100%",
            "value": 1.1507479861910241,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02063, 0.0198, 0.02], \"comparison_samples\": [1.256394640682095, 1.139240506329114, 1.145475372279496]}"
          },
          {
            "name": "igd/10K/basic/mid-gap2_100%",
            "value": 1.1293250141803743,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01991, 0.02014, 0.0198], \"comparison_samples\": [1.194361127774445, 1.142370958593307, 1.0921125206839493]}"
          },
          {
            "name": "igd/10K/basic/perfect-gap_100%",
            "value": 1.0844919786096257,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02179, 0.02028, 0.01998], \"comparison_samples\": [1.247995418098511, 1.0770047796070101, 1.0684491978609625]}"
          },
          {
            "name": "igd/10K/basic/perfect_100%",
            "value": 1.0629077353215284,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.02235, 0.02281, 0.02298], \"comparison_samples\": [1.0683556405353727, 1.0410771337288909, 1.070829450139795]}"
          },
          {
            "name": "igd/10K/basic/right-adjacent-gap_100%",
            "value": 1.0862445414847162,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.0199, 0.02017, 0.01966], \"comparison_samples\": [1.1909036505086776, 1.1009825327510918, 1.042970822281167]}"
          },
          {
            "name": "igd/10K/complex/mult_100bin",
            "value": 0.8746130030959752,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01149, 0.0113, 0.0113], \"comparison_samples\": [0.8893188854489165, 0.8841940532081377, 0.7718579234972677]}"
          },
          {
            "name": "intervaltree/10K/basic/3p-partial_100%",
            "value": 17.499485066941297,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.33859, 0.3409, 0.33984], \"comparison_samples\": [17.839304531085354, 17.554067971163747, 17.34762633996937]}"
          },
          {
            "name": "intervaltree/10K/basic/5p-partial_100%",
            "value": 17.346959632089934,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.34026, 0.33202, 0.33948], \"comparison_samples\": [18.323101777059776, 16.96576392437404, 16.98249124562281]}"
          },
          {
            "name": "intervaltree/10K/basic/contained_100%",
            "value": 14.274155164207519,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.2999, 0.30044, 0.29691], \"comparison_samples\": [14.665036674816626, 14.09854528390427, 14.131841980009519]}"
          },
          {
            "name": "intervaltree/10K/basic/enclosed_100%",
            "value": 16.80867232901207,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.38677, 0.37137, 0.37601], \"comparison_samples\": [18.243867924528303, 16.601251676352256, 16.24233261339093]}"
          },
          {
            "name": "intervaltree/10K/basic/left-adjacent-gap_100%",
            "value": 15.60888888888889,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.28096, 0.27671, 0.28857], \"comparison_samples\": [16.47859237536657, 15.364242087729039, 16.031666666666666]}"
          },
          {
            "name": "intervaltree/10K/basic/mid-gap1_100%",
            "value": 16.11680092059839,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.28233, 0.27966, 0.28011], \"comparison_samples\": [17.19427527405603, 16.090909090909093, 16.042955326460483]}"
          },
          {
            "name": "intervaltree/10K/basic/mid-gap2_100%",
            "value": 15.931934203062962,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.28088, 0.28432, 0.27776], \"comparison_samples\": [16.849430113977206, 16.127056154282474, 15.32046332046332]}"
          },
          {
            "name": "intervaltree/10K/basic/perfect-gap_100%",
            "value": 15.088770053475935,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.28216, 0.27026, 0.2847], \"comparison_samples\": [16.160366552119132, 14.35262878385555, 15.224598930481283]}"
          },
          {
            "name": "intervaltree/10K/basic/perfect_100%",
            "value": 15.525163094128612,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.33006, 0.34089, 0.33317], \"comparison_samples\": [15.777246653919693, 15.558649018712918, 15.525163094128612]}"
          },
          {
            "name": "intervaltree/10K/basic/right-adjacent-gap_100%",
            "value": 15.587882096069869,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.28557, 0.27233, 0.28621], \"comparison_samples\": [17.089766606822263, 14.865174672489085, 15.183554376657828]}"
          },
          {
            "name": "intervaltree/10K/complex/mult_100bin",
            "value": 31.229876160990717,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.41363, 0.4013, 0.40349], \"comparison_samples\": [32.01470588235294, 31.400625978090765, 27.560792349726775]}"
          },
          {
            "name": "polars_bio/10K/basic/3p-partial_100%",
            "value": 15.662718846549948,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.30462, 0.30417, 0.29958], \"comparison_samples\": [16.049525816649105, 15.662718846549948, 15.292496171516081]}"
          },
          {
            "name": "polars_bio/10K/basic/5p-partial_100%",
            "value": 15.465508431272355,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.30424, 0.30266, 0.29993], \"comparison_samples\": [16.3834141087776, 15.465508431272355, 15.004002001000499]}"
          },
          {
            "name": "polars_bio/10K/basic/contained_100%",
            "value": 14.316039980961445,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.29508, 0.34001, 0.30078], \"comparison_samples\": [14.429339853300734, 15.955419990614734, 14.316039980961445]}"
          },
          {
            "name": "polars_bio/10K/basic/enclosed_100%",
            "value": 13.776486365668305,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.30818, 0.33161, 0.29698], \"comparison_samples\": [14.53679245283019, 14.823871256146624, 12.828509719222463]}"
          },
          {
            "name": "polars_bio/10K/basic/left-adjacent-gap_100%",
            "value": 16.82277777777778,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.30281, 0.30535, 0.29056], \"comparison_samples\": [17.760117302052787, 16.95446973903387, 16.142222222222223]}"
          },
          {
            "name": "polars_bio/10K/basic/mid-gap1_100%",
            "value": 17.222669735327962,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.2986, 0.30889, 0.29933], \"comparison_samples\": [18.185140073081605, 17.772727272727273, 17.143757159221074]}"
          },
          {
            "name": "polars_bio/10K/basic/mid-gap2_100%",
            "value": 17.517867271695973,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.30884, 0.3155, 0.29793], \"comparison_samples\": [18.526694661067786, 17.895632444696542, 16.432984004412575]}"
          },
          {
            "name": "polars_bio/10K/basic/perfect-gap_100%",
            "value": 16.216042780748662,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.30324, 0.31529, 0.30249], \"comparison_samples\": [17.36769759450172, 16.744025491237387, 16.175935828877]}"
          },
          {
            "name": "polars_bio/10K/basic/perfect_100%",
            "value": 14.433830382106246,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.32527, 0.30874, 0.30975], \"comparison_samples\": [15.548279158699808, 14.091282519397536, 14.433830382106246]}"
          },
          {
            "name": "polars_bio/10K/basic/right-adjacent-gap_100%",
            "value": 16.203602620087338,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.29685, 0.31634, 0.29308], \"comparison_samples\": [17.764811490125673, 17.267467248908297, 15.548010610079578]}"
          },
          {
            "name": "polars_bio/10K/complex/mult_100bin",
            "value": 23.73606811145511,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.31427, 0.29798, 0.30667], \"comparison_samples\": [24.324303405572756, 23.31611893583725, 20.9474043715847]}"
          },
          {
            "name": "pygenogrove/10K/basic/3p-partial_100%",
            "value": 9.732234809474768,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.19101, 0.18563, 0.189], \"comparison_samples\": [10.063751317175976, 9.55870236869207, 9.647779479326188]}"
          },
          {
            "name": "pygenogrove/10K/basic/5p-partial_100%",
            "value": 9.50536535513541,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.19051, 0.18602, 0.18572], \"comparison_samples\": [10.259019924609586, 9.50536535513541, 9.29064532266133]}"
          },
          {
            "name": "pygenogrove/10K/basic/contained_100%",
            "value": 8.989528795811518,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.18623, 0.18917, 0.18887], \"comparison_samples\": [9.106601466992666, 8.877053026748007, 8.989528795811518]}"
          },
          {
            "name": "pygenogrove/10K/basic/enclosed_100%",
            "value": 8.416182387125614,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.18827, 0.19096, 0.18689], \"comparison_samples\": [8.88066037735849, 8.536432722396066, 8.073002159827213]}"
          },
          {
            "name": "pygenogrove/10K/basic/left-adjacent-gap_100%",
            "value": 8.87388888888889,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.16026, 0.15973, 0.15907], \"comparison_samples\": [9.399413489736071, 8.868961687951138, 8.837222222222222]}"
          },
          {
            "name": "pygenogrove/10K/basic/mid-gap1_100%",
            "value": 9.283659378596088,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.16063, 0.16159, 0.16135], \"comparison_samples\": [9.78258221680877, 9.29746835443038, 9.241122565864833]}"
          },
          {
            "name": "pygenogrove/10K/basic/mid-gap2_100%",
            "value": 9.247305728871243,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.16303, 0.16574, 0.15985], \"comparison_samples\": [9.779844031193761, 9.401020986954055, 8.816878102592387]}"
          },
          {
            "name": "pygenogrove/10K/basic/perfect-gap_100%",
            "value": 8.548663101604278,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.15792, 0.15986, 0.16279], \"comparison_samples\": [9.0446735395189, 8.489644184811471, 8.705347593582886]}"
          },
          {
            "name": "pygenogrove/10K/basic/perfect_100%",
            "value": 8.983224603914259,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.22079, 0.19278, 0.19058], \"comparison_samples\": [10.554015296367112, 8.798722044728436, 8.880708294501398]}"
          },
          {
            "name": "pygenogrove/10K/basic/right-adjacent-gap_100%",
            "value": 8.678493449781659,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.15899, 0.1597, 0.15845], \"comparison_samples\": [9.514661879114303, 8.717248908296943, 8.40583554376658]}"
          },
          {
            "name": "pygenogrove/10K/complex/mult_100bin",
            "value": 9.798761609907121,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.1266, 0.12995, 0.12608], \"comparison_samples\": [9.798761609907121, 10.168231611893585, 8.612021857923496]}"
          },
          {
            "name": "tabix/10K/basic/3p-partial_100%",
            "value": 11.23480947476828,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21818, 0.21279, 0.22556], \"comparison_samples\": [11.495258166491043, 10.957260556127704, 11.514037774374682]}"
          },
          {
            "name": "tabix/10K/basic/5p-partial_100%",
            "value": 11.22432294328053,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21966, 0.22327, 0.2178], \"comparison_samples\": [11.82875605815832, 11.408788962698006, 10.89544772386193]}"
          },
          {
            "name": "tabix/10K/basic/contained_100%",
            "value": 10.575440266539742,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.22219, 0.22207, 0.22824], \"comparison_samples\": [10.865036674816626, 10.42092914124824, 10.86339838172299]}"
          },
          {
            "name": "tabix/10K/basic/enclosed_100%",
            "value": 9.694233348234242,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21434, 0.21686, 0.2187], \"comparison_samples\": [10.110377358490567, 9.694233348234242, 9.44708423326134]}"
          },
          {
            "name": "tabix/10K/basic/left-adjacent-gap_100%",
            "value": 12.03777777777778,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21502, 0.21668, 0.22229], \"comparison_samples\": [12.611143695014663, 12.031093836757357, 12.349444444444444]}"
          },
          {
            "name": "tabix/10K/basic/mid-gap1_100%",
            "value": 12.451668584579977,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21641, 0.21242, 0.2194], \"comparison_samples\": [13.179658952496954, 12.222094361334868, 12.565864833906073]}"
          },
          {
            "name": "tabix/10K/basic/mid-gap2_100%",
            "value": 12.145207033465685,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21325, 0.21594, 0.21412], \"comparison_samples\": [12.79244151169766, 12.248440158820193, 11.810259238830668]}"
          },
          {
            "name": "tabix/10K/basic/perfect-gap_100%",
            "value": 11.444385026737967,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21401, 0.2098, 0.21619], \"comparison_samples\": [12.257159221076748, 11.141795007966012, 11.56096256684492]}"
          },
          {
            "name": "tabix/10K/basic/perfect_100%",
            "value": 10.225535880708295,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21944, 0.21921, 0.22017], \"comparison_samples\": [10.489483747609942, 10.005020538566864, 10.25955265610438]}"
          },
          {
            "name": "tabix/10K/basic/right-adjacent-gap_100%",
            "value": 11.832423580786026,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.21595, 0.21874, 0.21677], \"comparison_samples\": [12.923399162178336, 11.93995633187773, 11.49973474801061]}"
          },
          {
            "name": "tabix/10K/complex/mult_100bin",
            "value": 0.9102167182662539,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01165, 0.01176, 0.01194], \"comparison_samples\": [0.9017027863777091, 0.92018779342723, 0.8155737704918032]}"
          },
          {
            "name": "ucsc/10K/basic/3p-partial_100%",
            "value": 0.9335736354273945,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01813, 0.01806, 0.01846], \"comparison_samples\": [0.9552160168598525, 0.9299691040164778, 0.9423175089331292]}"
          },
          {
            "name": "ucsc/10K/basic/5p-partial_100%",
            "value": 0.8942258559018907,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01839, 0.01746, 0.0175], \"comparison_samples\": [0.9903069466882068, 0.8921819110884006, 0.8754377188594298]}"
          },
          {
            "name": "ucsc/10K/basic/contained_100%",
            "value": 0.8810090433127081,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01857, 0.01851, 0.01823], \"comparison_samples\": [0.9080684596577018, 0.8686062881276396, 0.8676820561637315]}"
          },
          {
            "name": "ucsc/10K/basic/enclosed_100%",
            "value": 0.8395172105498435,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01925, 0.01878, 0.01877], \"comparison_samples\": [0.9080188679245282, 0.8395172105498435, 0.8107991360691144]}"
          },
          {
            "name": "ucsc/10K/basic/left-adjacent-gap_100%",
            "value": 0.8733333333333335,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01572, 0.01641, 0.01534], \"comparison_samples\": [0.9219941348973608, 0.9111604664075513, 0.8522222222222222]}"
          },
          {
            "name": "ucsc/10K/basic/mid-gap1_100%",
            "value": 0.8883774453394707,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01544, 0.01611, 0.01507], \"comparison_samples\": [0.9403166869671132, 0.92692750287687, 0.8631156930126003]}"
          },
          {
            "name": "ucsc/10K/basic/mid-gap2_100%",
            "value": 0.8638684061259218,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01523, 0.0157, 0.01499], \"comparison_samples\": [0.913617276544691, 0.890527509926262, 0.8268063982349696]}"
          },
          {
            "name": "ucsc/10K/basic/perfect-gap_100%",
            "value": 0.8358288770053476,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01563, 0.01581, 0.01526], \"comparison_samples\": [0.8951890034364263, 0.8396176314391929, 0.816042780748663]}"
          },
          {
            "name": "ucsc/10K/basic/perfect_100%",
            "value": 0.8960857409133272,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01937, 0.01923, 0.0181], \"comparison_samples\": [0.925908221797323, 0.8776814240073026, 0.8434296365330849]}"
          },
          {
            "name": "ucsc/10K/basic/right-adjacent-gap_100%",
            "value": 0.8373362445414847,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01534, 0.01594, 0.01519], \"comparison_samples\": [0.9180131657690006, 0.8700873362445415, 0.8058355437665783]}"
          },
          {
            "name": "ucsc/10K/complex/mult_100bin",
            "value": 1.1222910216718267,
            "unit": "ratio to bedtools",
            "extra": "{\"samples_seconds\": [0.01472, 0.0145, 0.01416], \"comparison_samples\": [1.1393188854489165, 1.1345852895148671, 0.9672131147540984]}"
          },
          {
            "name": "ailist/all/end-to-end/wall",
            "value": 0.8740415180000127,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8541142600000171, 0.8740415180000127, 0.887100074999978], \"comparison_samples\": [0.8541142600000171, 0.8740415180000127, 0.887100074999978]}"
          },
          {
            "name": "awk/all/end-to-end/wall",
            "value": 15.825731768000026,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [15.743699942000006, 15.923285420000013, 15.825731768000026], \"comparison_samples\": [15.743699942000006, 15.923285420000013, 15.825731768000026]}"
          },
          {
            "name": "bedmaps/all/end-to-end/wall",
            "value": 1.0147909780000077,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.9819043569999906, 1.0147909780000077, 1.0274569000000042], \"comparison_samples\": [0.9819043569999906, 1.0147909780000077, 1.0274569000000042]}"
          },
          {
            "name": "bedops/all/end-to-end/wall",
            "value": 0.9734507099999519,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.9463154400000064, 0.9734507099999519, 0.9812803199999962], \"comparison_samples\": [0.9463154400000064, 0.9734507099999519, 0.9812803199999962]}"
          },
          {
            "name": "bedtk/all/end-to-end/wall",
            "value": 0.788256959000023,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.791321757999981, 0.788256959000023, 0.7850861599999917], \"comparison_samples\": [0.791321757999981, 0.788256959000023, 0.7850861599999917]}"
          },
          {
            "name": "bedtk_sorted/all/end-to-end/wall",
            "value": 0.9275380720000044,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.9130112959999792, 0.9275380720000044, 0.9300929710000219], \"comparison_samples\": [0.9130112959999792, 0.9275380720000044, 0.9300929710000219]}"
          },
          {
            "name": "bedtools/all/end-to-end/wall",
            "value": 0.9071141059999945,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8854388819999883, 0.9071141059999945, 0.925813536000021], \"comparison_samples\": [0.8854388819999883, 0.9071141059999945, 0.925813536000021]}"
          },
          {
            "name": "bedtools_sorted/all/end-to-end/wall",
            "value": 0.9663971370000013,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.9425685419999752, 0.9764641330000359, 0.9663971370000013], \"comparison_samples\": [0.9425685419999752, 0.9764641330000359, 0.9663971370000013]}"
          },
          {
            "name": "genogrove/all/end-to-end/wall",
            "value": 0.8818192389999808,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.9579332049999607, 0.8545564130000116, 0.8818192389999808], \"comparison_samples\": [0.9579332049999607, 0.8545564130000116, 0.8818192389999808]}"
          },
          {
            "name": "gia/all/end-to-end/wall",
            "value": 0.851966066999978,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.8868577060000007, 0.8459057119999898, 0.851966066999978], \"comparison_samples\": [0.8868577060000007, 0.8459057119999898, 0.851966066999978]}"
          },
          {
            "name": "giggle/all/end-to-end/wall",
            "value": 8.536268813999982,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [8.675736907000015, 8.536268813999982, 8.509774094000022], \"comparison_samples\": [8.675736907000015, 8.536268813999982, 8.509774094000022]}"
          },
          {
            "name": "granges/all/end-to-end/wall",
            "value": 0.8540014779999865,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.9694893889999889, 0.8540014779999865, 0.8170045720000303], \"comparison_samples\": [0.9694893889999889, 0.8540014779999865, 0.8170045720000303]}"
          },
          {
            "name": "igd/all/end-to-end/wall",
            "value": 1.349937347999969,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [1.3235448179999878, 1.349937347999969, 1.3545840390000308], \"comparison_samples\": [1.3235448179999878, 1.349937347999969, 1.3545840390000308]}"
          },
          {
            "name": "intervaltree/all/end-to-end/wall",
            "value": 4.202208525000003,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [4.2207351730000084, 4.160620220999988, 4.202208525000003], \"comparison_samples\": [4.2207351730000084, 4.160620220999988, 4.202208525000003]}"
          },
          {
            "name": "polars_bio/all/end-to-end/wall",
            "value": 4.005485470999986,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [4.005485470999986, 4.0916222910000215, 3.93192072100004], \"comparison_samples\": [4.005485470999986, 4.0916222910000215, 3.93192072100004]}"
          },
          {
            "name": "pygenogrove/all/end-to-end/wall",
            "value": 2.5174202700000023,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [2.5554249809999874, 2.5174202700000023, 2.500390698999979], \"comparison_samples\": [2.5554249809999874, 2.5174202700000023, 2.500390698999979]}"
          },
          {
            "name": "segmeter/all/simulation/wall",
            "value": 1.4489666420000162,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [1.4489666420000162], \"comparison_samples\": [1.4489666420000162]}"
          },
          {
            "name": "tabix/all/end-to-end/wall",
            "value": 2.9384818749999795,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [2.901457327000003, 2.9385263080000072, 2.9384818749999795], \"comparison_samples\": [2.901457327000003, 2.9385263080000072, 2.9384818749999795]}"
          },
          {
            "name": "ucsc/all/end-to-end/wall",
            "value": 0.9054269160000104,
            "unit": "seconds",
            "extra": "{\"samples_seconds\": [0.9054269160000104, 0.9067831870000305, 0.8929392300000245], \"comparison_samples\": [0.9054269160000104, 0.9067831870000305, 0.8929392300000245]}"
          }
        ]
      }
    ]
  }
}