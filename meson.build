project('fmm2d', ['fortran', 'c'])

pkg = import('pkgconfig')

fmm2d_src = [
  # ./src/biharmonic
  'src/biharmonic/bh2dterms.f',
  'src/biharmonic/bhfmm2d.f',
  'src/biharmonic/bhfmm2dwrap.f',
  'src/biharmonic/bhkernels2d.f',
  'src/biharmonic/bhndiv2d.f',
  'src/biharmonic/bhrouts2d.f',

  # ./src/common
  'src/common/cdjseval2d.f',
  'src/common/cumsum.f',
  # 'src/common/dfft.f', # TODO: unused?
  'src/common/dfft_threadsafe.f',
  'src/common/dlaran.f',
  'src/common/fmmcommon2d.f',
  'src/common/hank103.f',
  'src/common/hkrand.f',
  'src/common/next235.f',
  'src/common/prini.f',
  'src/common/pts_tree2d.f',
  'src/common/tree_routs2d.f',

  # ./src/helmholtz
  'src/helmholtz_c.f03', # c wrappers
  'src/helmholtz/h2dcommon.f',
  'src/helmholtz/h2dterms.f',
  'src/helmholtz/helmkernels2d.f',
  'src/helmholtz/helmrouts2d.f',
  'src/helmholtz/hfmm2d.f',
  'src/helmholtz/hfmm2d_mps.f',
  'src/helmholtz/hfmm2d_mps_ndiv.f',
  'src/helmholtz/hfmm2d_ndiv.f',
  'src/helmholtz/hfmm2dwrap.f',
  'src/helmholtz/hfmm2dwrap_vec.f',
  'src/helmholtz/hndiv2d.f',
  'src/helmholtz/wideband2d.f',

  # ./src/laplace
  'src/laplace/cauchykernels2d.f',
  'src/laplace/cfmm2d.f',
  'src/laplace/cfmm2d_ndiv.f',
  'src/laplace/cfmm2dwrap.f',
  'src/laplace/cfmm2dwrap_vec.f',
  'src/laplace/l2dterms.f',
  'src/laplace/lapkernels2d.f',
  'src/laplace/laprouts2d.f',
  'src/laplace/lfmm2d.f',
  'src/laplace/lfmm2d_ndiv.f',
  'src/laplace/lfmm2dwrap.f',
  'src/laplace/lfmm2dwrap_vec.f',
  'src/laplace/lndiv2d.f',
  'src/laplace/rfmm2d.f',
  'src/laplace/rfmm2d_ndiv.f',
  'src/laplace/rfmm2dwrap.f',
  'src/laplace/rfmm2dwrap_vec.f',
  'src/laplace/rlapkernels2d.f',

  # TODO: there's some code generation using happening with the
  # modified biharmonic kernel (commented out lines below). Look into
  # this later.

  # ./src/modified-biharmonic
  # 'src/modified-biharmonic/jinjaroot.yaml',
  # 'src/modified-biharmonic/jinjaroot.yaml.py',
  # 'src/modified-biharmonic/mbhfmm2d.f',
  # 'src/modified-biharmonic/mbhgreen2d.f',
  # 'src/modified-biharmonic/mbhkernels2d.f',
  # 'src/modified-biharmonic/mbhkernels2d.f.j2',
  # 'src/modified-biharmonic/mbhrouts2d.f',
  # 'src/modified-biharmonic/mbhrouts2d.f.j2',

  # ./src/stokes
  'src/stokes/stfmm2d.f',
  'src/stokes/stokkernels2d.f',
]

# TODO: set these the right way using meson instead of doing it manually...
fmm2d_fortran_args = ['-fPIC', '-O3', '-march=native', '-funroll-loops', '-std=legacy', '-w']
fmm2d_c_args = ['-std=c99'] + fmm2d_fortran_args

fmm2d_lib = library(
  'fmm2d',
  fmm2d_src,
  fortran_args : fmm2d_fortran_args,
  install : true
)

install_subdir('include/fmm2d', install_dir : 'include', install_tag : 'devel')

pkg.generate(fmm2d_lib)
