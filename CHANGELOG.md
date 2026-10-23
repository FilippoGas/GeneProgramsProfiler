# Changelog

## [1.2.0](https://github.com/FilippoGas/GeneProgramsProfiler/compare/v1.1.0...v1.2.0) (2026-10-23)


### Features

* 6 : add report captions ([02d0e2f](https://github.com/FilippoGas/GeneProgramsProfiler/commit/02d0e2f2187b72739f1367769fc09731d9ca8055))
* 7 : Enabled rules benchmarking ([bdae44c](https://github.com/FilippoGas/GeneProgramsProfiler/commit/bdae44c138a2543fd18cfdc3a780e31b82214556))
* add comparative table caption and fix category naming in make_comp_table rule ([8d871b7](https://github.com/FilippoGas/GeneProgramsProfiler/commit/8d871b7c3f01d8fe18ae26ecdacbce9d6de0f4f8))
* add logging to cNMF workflow rules for better traceability ([d120233](https://github.com/FilippoGas/GeneProgramsProfiler/commit/d120233e758bc242ebf3864ebf31127de140a3fa))
* add r-bpcells dependency to DE_analysis and preprocess environments ([19c8a80](https://github.com/FilippoGas/GeneProgramsProfiler/commit/19c8a803e3e29443e6f5581bab66a2d03d230ff2))
* add vertical upset intersection table and update report captions ([f958845](https://github.com/FilippoGas/GeneProgramsProfiler/commit/f958845c61016e0f3d0f2884cded9173f01f5f5e))
* enhance k selection process in extract_best_k.R with validation and delta computation ([290799c](https://github.com/FilippoGas/GeneProgramsProfiler/commit/290799ce18178575ce5a87169394fb1dc2d60dda))
* enhance memory management and export functionality in annotate_and_save.R ([c569ef1](https://github.com/FilippoGas/GeneProgramsProfiler/commit/c569ef10d36a7dccd1798be42b6ff8aaf7d8d62e))
* update gene list extraction and add BPCells library in run_DE_analysis.R ([97fbd16](https://github.com/FilippoGas/GeneProgramsProfiler/commit/97fbd16d02346b757f196a78986ad45a0fdbc60d))


### Bug Fixes

* cell_id column was missing from metadata ([e50e053](https://github.com/FilippoGas/GeneProgramsProfiler/commit/e50e053947b0e18263d8771b9b9f5b4c0cd42afb))
* Correct metadata path ([92569bf](https://github.com/FilippoGas/GeneProgramsProfiler/commit/92569bfc3478b5eb36c0efbd791f942629363ca6))
* correct report syntax for horizontal and vertical upset plots ([96c3a31](https://github.com/FilippoGas/GeneProgramsProfiler/commit/96c3a31408852bdb78c6dbcb31bc94ee79b2b3de))
* ensure r-bpcells dependency is correctly formatted in preprocess.yaml ([76f59c5](https://github.com/FilippoGas/GeneProgramsProfiler/commit/76f59c571de0bae3d1e174ee16e3aba13e4ab0b4))
* fixed metadata columns before save ([8dfc37a](https://github.com/FilippoGas/GeneProgramsProfiler/commit/8dfc37acbd2b167c14cd334c961dec9e81259f4c))
* update benchmark file naming for cNMF factorization worker ([c5821fe](https://github.com/FilippoGas/GeneProgramsProfiler/commit/c5821fe2fe14962dc4704086b21df882b69afdfd))

## [1.1.0](https://github.com/FilippoGas/GeneProgramsProfiler/compare/v1.0.0...v1.1.0) (2026-09-16)


### Features

* add release configuration and manifest files ([76ee474](https://github.com/FilippoGas/GeneProgramsProfiler/commit/76ee47495e8a3379c7a5af90b9ef084d17fd7029))


### Bug Fixes

* update job status reporting in PBS and Slurm scripts ([f022940](https://github.com/FilippoGas/GeneProgramsProfiler/commit/f022940c0b0086b1c45dcfb64356a3a06184bf36))

## [1.0.0](https://github.com/FilippoGas/GeneProgramsProfiler/releases/tag/v1.0.0)

Initial release of GeneProgramsProfiler.
