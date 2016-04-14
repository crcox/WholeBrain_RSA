require 'rake'
require 'rake/clean'

HOME = ENV['HOME']
PERMDIR = "../../../permtest_rows/sample/solutionmaps"
SHARED_ATLAS = "#{HOME}/MRI/Manchester/data/CommonBrains/MNI_EPI_funcRes.nii"
SPEC_BOTH = "#{HOME}/suma_TT_N27/TT_N27_both.spec"
SURFACE_VOLUME = "./TT_N27_SurfVol.nii"

TXT = Rake::FileList["txt/??.mni"]
TXT_IGRAPH = Rake::FileList["txt/??_igraph.mni"]
TXT_ABS = Rake::FileList["txt/??_abs.mni"]
TXT_MASK = Rake::FileList["txt/??_mask.mni"]

AFNI_RAW = TXT.pathmap("afni/%n_raw+tlrc.HEAD")
AFNI = TXT.pathmap("afni/%n+tlrc.HEAD")
AFNIBLUR = TXT.pathmap("afni/%n.b4+tlrc.HEAD")
AFNI_IGRAPH = TXT_IGRAPH.pathmap("afni/%n+tlrc.HEAD")
AFNI_ABS= TXT_ABS.pathmap("afni/%n+tlrc.HEAD")
AFNI_MASK = TXT_MASK.pathmap("afni/%n+tlrc.HEAD")

RANK = AFNI.pathmap("rank/%f")
RANKBLUR = RANK.sub("+tlrc",".b4+tlrc")

ZSCORE = AFNI.pathmap("zscore/%f")
ZSCORE_IGRAPH = AFNI_IGRAPH.pathmap("zscore/%f")
ZSCORE_ABS = AFNI_ABS.pathmap("zscore/%f")
ZSCOREBLUR = ZSCORE.sub("+tlrc",".b4+tlrc")
ZSCOREBLUR_IGRAPH = ZSCORE_IGRAPH.sub("+tlrc",".b4+tlrc")
ZSCOREBLUR_ABS = ZSCORE_ABS.sub("+tlrc",".b4+tlrc")

TTEST = "ttest_nodestrength+tlrc.HEAD"
TTEST_IGRAPH = "ttest_nodestrength_igraph+tlrc.HEAD"
TTEST_ABS = "ttest_nodestrength_abs+tlrc.HEAD"
TTESTBLUR = "ttest_nodestrength.b4+tlrc.HEAD"
TTESTBLUR_IGRAPH = "ttest_nodestrength_igraph.b4+tlrc.HEAD"
TTESTBLUR_ABS = "ttest_nodestrength_abs.b4+tlrc.HEAD"

NONPARAMETRIC = "nonparametric_nodestrength+tlrc.HEAD"
NONPARAMETRIC_MEAN = "nonparametric_meanrank+tlrc.HEAD"

AVG = "mean_nodestrength+tlrc.HEAD"
AVG_RANK = "mean_rank+tlrc.HEAD"
OVERLAP = "overlap+tlrc.HEAD"
PERM_MASK = "permtest_mask+tlrc.HEAD"
GROUP_MASK = "group_mask+tlrc.HEAD"
PERM = "#{PERMDIR}/permutations+tlrc.HEAD"
PERMMEAN = AFNI.pathmap("#{PERMDIR}/avg/subj/%f")
PERMSD = AFNI.pathmap("#{PERMDIR}/sd/subj/%f")

PNG_TTEST = [TTEST.sub('+tlrc.HEAD','_p05.png'),TTEST.sub('+tlrc.HEAD','_p01.png')]
PNG_TTESTBLUR = [TTESTBLUR.sub('+tlrc.HEAD','_p05.png'),TTESTBLUR.sub('+tlrc.HEAD','_p01.png')]
PNG_NONPARAMETRIC = [NONPARAMETRIC.sub('+tlrc.HEAD','_p05.png'),NONPARAMETRIC.sub('+tlrc.HEAD','_p01.png'),NONPARAMETRIC.sub('+tlrc.HEAD','_p001.png')]
PNG_NONPARAMETRIC_MEAN = [NONPARAMETRIC_MEAN.sub('+tlrc.HEAD','_p05.png'),NONPARAMETRIC_MEAN.sub('+tlrc.HEAD','_p01.png'),NONPARAMETRIC_MEAN.sub('+tlrc.HEAD','_p001.png')]

task :start_afni do
  system("afni -niml -yesplugouts &")
  sh("plugout_drive -com 'SWITCH_ANATOMY #{SURFACE_VOLUME}' -quit")
  sh("plugout_drive -com 'SET_THRESHNEW 0' -quit")
  sh("plugout_drive -com 'SET_PBAR_SIGN -' -quit")
  sh("plugout_drive -com 'SEE_OVERLAY +' -quit")

  system("suma -niml -spec #{SPEC_BOTH} -sv #{SURFACE_VOLUME} &")
  sh("DriveSuma -com  viewer_cont -key:d 't'")           # talk to afni
  sh("DriveSuma -com  viewer_cont -key:r2 '.'")          # Select the inflated surfaces
  sh("DriveSuma -com  viewer_cont -key 'F3'")            # toggle the crosshair (off)
  sh("DriveSuma -com  viewer_cont -key 'F6'")            # toggle background color (to white from black)
  sh("DriveSuma -com  viewer_cont -key 'F9'")            # toggle the label at the crosshair (off)
  sh("DriveSuma -com  viewer_cont -viewer_size 700 600") # adjust viewer size (which effects figure size)
end

task :stop_afni do
  sh("DriveSuma -com kill_suma")
  sh("plugout_drive -com 'QUIT' -quit")
end

PNG_TTEST.each do |target|
  file target => TTEST do
    prefix = TTEST.sub('+tlrc.HEAD','')
    sh("plugout_drive -com 'RESCAN_THIS' -quit")
    sh("plugout_drive -com 'SET_FUNCTION #{TTEST}' -quit")
    if target.include? "p05" then
      sh("plugout_drive -com 'SET_THRESHNEW A 0.05 *p' -quit")
    elsif target.include? "p01" then
      sh("plugout_drive -com 'SET_THRESHNEW A 0.01 *p' -quit")
    end
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.left.ppm' -key:d 'ctrl+left' -key 'ctrl+r'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.right_medial.ppm' -key:d '[' -key 'ctrl+r' -key:d '['")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.bottom.ppm' -key:d 'ctrl+down' -key 'ctrl+r'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.right.ppm' -key:d 'ctrl+right' -key 'ctrl+r'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.left_medial.ppm' -key:d ']' -key 'ctrl+r' -key:d ']'")
    sh("DriveSuma -com viewer_cont -key:r17 'right'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.back.ppm' -key:d 'right' -key 'ctrl+r'")

    sh("imcat -prefix #{target.sub('.png','')} -nx 3 -ny 2 #{prefix}.left.*.ppm #{prefix}.back.*.ppm #{prefix}.right.*.ppm #{prefix}.left_medial.*.ppm #{prefix}.bottom.*.ppm #{prefix}.right_medial.*.ppm")
    sh("mogrify -format png *.ppm")
    sh("rm *.ppm")
  end
  tmp_png = Rake::FileList["#{TTEST.sub('+tlrc.HEAD','')}.*.png"].exclude(target)
  CLEAN.push(tmp_png)
  CLOBBER.push(target)
end

PNG_TTESTBLUR.each do |target|
  file target => TTESTBLUR do
    prefix = TTESTBLUR.sub('+tlrc.HEAD','')
    sh("plugout_drive -com 'RESCAN_THIS' -quit")
    sh("plugout_drive -com 'SET_FUNCTION #{TTESTBLUR}' -quit")
    if target.include? "p05" then
      p="p05"
      sh("plugout_drive -com 'SET_THRESHNEW A 0.05 *p' -quit")
    elsif target.include? "p01" then
      p="p01"
      sh("plugout_drive -com 'SET_THRESHNEW A 0.01 *p' -quit")
    end
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.left.ppm' -key:d 'ctrl+left' -key 'ctrl+r'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.right_medial.ppm' -key:d '[' -key 'ctrl+r' -key:d '['")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.bottom.ppm' -key:d 'ctrl+down' -key 'ctrl+r'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.right.ppm' -key:d 'ctrl+right' -key 'ctrl+r'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.left_medial.ppm' -key:d ']' -key 'ctrl+r' -key:d ']'")
    sh("DriveSuma -com viewer_cont -key:r17 'right'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.back.ppm' -key:d 'right' -key 'ctrl+r'")

    sh("imcat -prefix #{target.sub('.png','')} -nx 3 -ny 2 #{prefix}.left.*.ppm #{prefix}.back.*.ppm #{prefix}.right.*.ppm #{prefix}.left_medial.*.ppm #{prefix}.bottom.*.ppm #{prefix}.right_medial.*.ppm")
    sh("mogrify -format png *.ppm")
    sh("rm *.ppm")
  end
  tmp_png = Rake::FileList["#{TTESTBLUR.sub('+tlrc.HEAD','')}.*.png"].exclude(target)
  CLEAN.push(tmp_png)
  CLOBBER.push(target)
end

PNG_NONPARAMETRIC.each do |target|
  file target => NONPARAMETRIC do
    prefix = NONPARAMETRIC.sub('+tlrc.HEAD','')
    sh("plugout_drive -com 'RESCAN_THIS' -quit")
    sh("plugout_drive -com 'SET_FUNCTION #{NONPARAMETRIC}' -quit")
    if target.include? "p05" then
      p="p05"
      sh("plugout_drive -com 'SET_THRESHNEW A 16 *' -quit")
    elsif target.include? "p01" then
      p="p01"
      sh("plugout_drive -com 'SET_THRESHNEW A 18 *' -quit")
    elsif target.include? "p001" then
      p="p001"
      sh("plugout_drive -com 'SET_THRESHNEW A 20 *' -quit")
    end
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.left.ppm' -key:d 'ctrl+left' -key 'ctrl+r'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.right_medial.ppm' -key:d '[' -key 'ctrl+r' -key:d '['")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.bottom.ppm' -key:d 'ctrl+down' -key 'ctrl+r'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.right.ppm' -key:d 'ctrl+right' -key 'ctrl+r'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.left_medial.ppm' -key:d ']' -key 'ctrl+r' -key:d ']'")
    sh("DriveSuma -com viewer_cont -key:r17 'right'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.back.ppm' -key:d 'right' -key 'ctrl+r'")

    sh("imcat -prefix #{target.sub('.png','')} -nx 3 -ny 2 #{prefix}.left.*.ppm #{prefix}.back.*.ppm #{prefix}.right.*.ppm #{prefix}.left_medial.*.ppm #{prefix}.bottom.*.ppm #{prefix}.right_medial.*.ppm")
    sh("mogrify -format png *.ppm")
    sh("rm *.ppm")
  end
  tmp_png = Rake::FileList["#{NONPARAMETRIC.sub('+tlrc.HEAD','')}.*.png"].exclude(target)
  CLEAN.push(tmp_png)
  CLOBBER.push(target)
end

PNG_NONPARAMETRIC_MEAN.each do |target|
  file target => NONPARAMETRIC_MEAN do
    prefix = NONPARAMETRIC_MEAN.sub('+tlrc.HEAD','')
    sh("plugout_drive -com 'RESCAN_THIS' -quit")
    sh("plugout_drive -com 'SET_FUNCTION #{NONPARAMETRIC_MEAN}' -quit")
    if target.include? "p05" then
      p="p05"
      sh("plugout_drive -com 'SET_THRESHNEW A 16 *' -quit")
    elsif target.include? "p01" then
      p="p01"
      sh("plugout_drive -com 'SET_THRESHNEW A 18 *' -quit")
    elsif target.include? "p001" then
      p="p001"
      sh("plugout_drive -com 'SET_THRESHNEW A 20 *' -quit")
    end
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.left.ppm' -key:d 'ctrl+left' -key 'ctrl+r'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.right_medial.ppm' -key:d '[' -key 'ctrl+r' -key:d '['")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.bottom.ppm' -key:d 'ctrl+down' -key 'ctrl+r'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.right.ppm' -key:d 'ctrl+right' -key 'ctrl+r'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.left_medial.ppm' -key:d ']' -key 'ctrl+r' -key:d ']'")
    sh("DriveSuma -com viewer_cont -key:r17 'right'")
    sh("DriveSuma -com viewer_cont -autorecord '#{prefix}.back.ppm' -key:d 'right' -key 'ctrl+r'")

    sh("imcat -prefix #{target.sub('.png','')} -nx 3 -ny 2 #{prefix}.left.*.ppm #{prefix}.back.*.ppm #{prefix}.right.*.ppm #{prefix}.left_medial.*.ppm #{prefix}.bottom.*.ppm #{prefix}.right_medial.*.ppm")
    sh("mogrify -format png *.ppm")
    sh("rm *.ppm")
  end
  tmp_png = Rake::FileList["#{NONPARAMETRIC_MEAN.sub('+tlrc.HEAD','')}.*.png"].exclude(target)
  CLEAN.push(tmp_png)
  CLOBBER.push(target)
end

AFNI_RAW.zip(TXT).each do |target,source|
  file target => [source] do
    sh("3dUndump -master #{SHARED_ATLAS} -xyz -datum float -prefix #{target.sub("+tlrc.HEAD","")} #{source}")
  end
  CLOBBER.push(target)
  CLOBBER.push(target.sub(".HEAD",".BRIK"))
  CLOBBER.push(target.sub(".HEAD",".BRIK.gz"))
end

AFNI.zip(AFNI_RAW).each do |target,source|
  file target => [source] do
    sh("3dcalc -a #{source} -expr 'a*1000' -prefix #{target.sub("+tlrc.HEAD","")}")
  end
  CLOBBER.push(target)
  CLOBBER.push(target.sub(".HEAD",".BRIK"))
  CLOBBER.push(target.sub(".HEAD",".BRIK.gz"))
end

AFNIBLUR.zip(AFNI).each do |target,source|
  file target => [source] do
    puts "#{target} => #{source}"
    sh("3dmerge -1blur_fwhm 4 -prefix #{target.sub('+tlrc.HEAD','')} #{source}")
  end
  CLOBBER.push(target)
  CLOBBER.push(target.sub(".HEAD",".BRIK"))
  CLOBBER.push(target.sub(".HEAD",".BRIK.gz"))
end

AFNI_MASK.zip(TXT_MASK).each do |target,source|
  file target => [source] do
    sh("3dUndump -master #{SHARED_ATLAS} -xyz -datum float -prefix #{target.sub("+tlrc.HEAD","")} #{source}")
  end
  CLOBBER.push(target)
  CLOBBER.push(target.sub(".HEAD",".BRIK"))
  CLOBBER.push(target.sub(".HEAD",".BRIK.gz"))
end

AFNI_ABS.zip(TXT_ABS).each do |target,source|
  file target => [source] do
    sh("3dUndump -master #{SHARED_ATLAS} -xyz -datum float -prefix #{target.sub("+tlrc.HEAD","")} #{source}")
  end
  CLOBBER.push(target)
  CLOBBER.push(target.sub(".HEAD",".BRIK"))
  CLOBBER.push(target.sub(".HEAD",".BRIK.gz"))
end

AFNI_IGRAPH.zip(TXT_IGRAPH).each do |target,source|
  file target => [source] do
    sh("3dUndump -master #{SHARED_ATLAS} -xyz -datum float -prefix #{target.sub("+tlrc.HEAD","")} #{source}")
  end
  CLOBBER.push(target)
  CLOBBER.push(target.sub(".HEAD",".BRIK"))
  CLOBBER.push(target.sub(".HEAD",".BRIK.gz"))
end

file AVG => AFNI do
  sh("3dmerge -1blur_fwhm 4 -gmean -prefix #{AVG.sub('+tlrc.HEAD','')} #{AFNI.join(' ')}")
end
CLOBBER.push(AVG)
CLOBBER.push(AVG.sub(".HEAD",".BRIK"))
CLOBBER.push(AVG.sub(".HEAD",".BRIK.gz"))

file AVG_RANK => RANK do
  sh("3dmerge -1blur_fwhm 4 -gmean -prefix #{AVG_RANK.sub('+tlrc.HEAD','')} #{RANK.join(' ')}")
end
CLOBBER.push(AVG_RANK)
CLOBBER.push(AVG_RANK.sub(".HEAD",".BRIK"))
CLOBBER.push(AVG_RANK.sub(".HEAD",".BRIK.gz"))

RANK.zip(AFNI).each do |target,source|
  file target => source do
    puts "#{target} => #{source}"
    basename = File.basename(target,'+tlrc.HEAD')
    pattern = "#{PERMDIR}/afni/???_#{basename}+tlrc.HEAD"
    permlist = Rake::FileList[pattern]
    rankList = []
    permlist.each_with_index do |permutation,iPerm|
      rankPrefix = target.sub('+tlrc.HEAD',"_#{iPerm}")
      rankList.push(rankPrefix)
      sh("3dcalc -a #{source} -b #{permutation} -expr 'step(a-b)' -prefix #{rankPrefix}")
    end
    sh("3dmerge -gcount -prefix #{target.sub('+tlrc.HEAD','')} #{rankList.collect {|x| x+'+tlrc.HEAD'}.join(' ')}")
    rankList.each do |prefix|
      rm "#{prefix}+tlrc.HEAD"
      rm "#{prefix}+tlrc.BRIK.gz"
    end
  end
  CLOBBER.push(target)
  CLOBBER.push(target.sub(".HEAD",".BRIK"))
  CLOBBER.push(target.sub(".HEAD",".BRIK.gz"))
end

RANKBLUR.zip(RANK).each do |target,source|
  file target => source do
    sh("3dmerge -1blur_fwhm 4 -prefix #{target.sub('+tlrc.HEAD','')} #{source}")
  end
  CLOBBER.push(target)
  CLOBBER.push(target.sub(".HEAD",".BRIK"))
  CLOBBER.push(target.sub(".HEAD",".BRIK.gz"))
end

file NONPARAMETRIC => RANK+[AVG] do
  sh("3dmerge -1clip 50 -1blur_fwhm 4 -gcount -prefix rm_rankcount #{RANK.join(' ')}")

  sh("3dbucket -fbuc -prefix #{NONPARAMETRIC.sub('+tlrc.HEAD','')} #{AVG} rm_rankcount+tlrc.HEAD")
  sh("3drefit -fith #{NONPARAMETRIC}")

  CLEAN.push("rm_rankcount+tlrc.HEAD")
  CLEAN.push("rm_rankcount+tlrc.BRIK")
  CLEAN.push("rm_rankcount+tlrc.BRIK.gz")
end
CLOBBER.push(NONPARAMETRIC)
CLOBBER.push(NONPARAMETRIC.sub(".HEAD",".BRIK"))
CLOBBER.push(NONPARAMETRIC.sub(".HEAD",".BRIK.gz"))

file NONPARAMETRIC_MEAN => RANK+[AVG_RANK] do
  sh("3dmerge -1clip 50 -1blur_fwhm 4 -gcount -prefix rm_rankcount #{RANK.join(' ')}")

  sh("3dbucket -fbuc -prefix #{NONPARAMETRIC_MEAN.sub('+tlrc.HEAD','')} #{AVG_RANK} rm_rankcount+tlrc.HEAD")
  sh("3drefit -fith #{NONPARAMETRIC_MEAN}")

  CLEAN.push("rm_rankcount+tlrc.HEAD")
  CLEAN.push("rm_rankcount+tlrc.BRIK")
  CLEAN.push("rm_rankcount+tlrc.BRIK.gz")
end
CLOBBER.push(NONPARAMETRIC_MEAN)
CLOBBER.push(NONPARAMETRIC_MEAN.sub(".HEAD",".BRIK"))
CLOBBER.push(NONPARAMETRIC_MEAN.sub(".HEAD",".BRIK.gz"))

ZSCORE.zip(AFNI,PERMMEAN,PERMSD) do |z,a,pm,psd|
  file z => [a,pm,psd] do
    sh("3dcalc -a #{a} -b #{pm} -c #{psd} -expr 'min((a-b)/c,5)*notzero(a)' -prefix #{z.sub("+tlrc.HEAD","")}")
  end
  CLOBBER.push(z)
  CLOBBER.push(z.sub(".HEAD",".BRIK"))
  CLOBBER.push(z.sub(".HEAD",".BRIK.gz"))
end

ZSCORE_ABS.zip(AFNI_ABS,PERMMEAN,PERMSD) do |z,a,pm,psd|
  file z => [a,pm,psd] do
    sh("3dcalc -a #{a} -b #{pm} -c #{psd} -expr 'min((a-b)/c,5)*notzero(a)' -prefix #{z.sub("+tlrc.HEAD","")}")
  end
  CLOBBER.push(z)
  CLOBBER.push(z.sub(".HEAD",".BRIK"))
  CLOBBER.push(z.sub(".HEAD",".BRIK.gz"))
end

ZSCORE_IGRAPH.zip(AFNI_IGRAPH,PERMMEAN,PERMSD) do |z,a,pm,psd|
  file z => [a,pm,psd] do
    sh("3dcalc -a #{a} -b #{pm} -c #{psd} -expr 'min((a-b)/c,5)*notzero(a)' -prefix #{z.sub("+tlrc.HEAD","")}")
  end
  CLOBBER.push(z)
  CLOBBER.push(z.sub(".HEAD",".BRIK"))
  CLOBBER.push(z.sub(".HEAD",".BRIK.gz"))
end

ZSCOREBLUR.zip(ZSCORE) do |target,source|
  file target => source do
    sh("3dmerge -1blur_fwhm 4.0 -prefix #{target.sub("+tlrc.HEAD","")} #{source}")
  end
  CLOBBER.push(target)
  CLOBBER.push(target.sub(".HEAD",".BRIK"))
  CLOBBER.push(target.sub(".HEAD",".BRIK.gz"))
end

ZSCOREBLUR_ABS.zip(ZSCORE_ABS) do |target,source|
  file target => source do
    sh("3dmerge -1blur_fwhm 4.0 -prefix #{target.sub("+tlrc.HEAD","")} #{source}")
  end
  CLOBBER.push(target)
  CLOBBER.push(target.sub(".HEAD",".BRIK"))
  CLOBBER.push(target.sub(".HEAD",".BRIK.gz"))
end

ZSCOREBLUR_IGRAPH.zip(ZSCORE_IGRAPH) do |target,source|
  file target => source do
    sh("3dmerge -1blur_fwhm 4.0 -prefix #{target.sub("+tlrc.HEAD","")} #{source}")
  end
  CLOBBER.push(target)
  CLOBBER.push(target.sub(".HEAD",".BRIK"))
  CLOBBER.push(target.sub(".HEAD",".BRIK.gz"))
end

file GROUP_MASK => AFNI_MASK do
  sh("3dmerge -gcount -prefix #{GROUP_MASK.sub("+tlrc.HEAD","")} #{AFNI_MASK.join(" ")}")
end
CLOBBER.push(GROUP_MASK)
CLOBBER.push(GROUP_MASK.sub(".HEAD",".BRIK"))
CLOBBER.push(GROUP_MASK.sub(".HEAD",".BRIK.gz"))

file TTEST=> ZSCORE+[GROUP_MASK] do
  sh("3dttest++ -mask #{GROUP_MASK} -setA #{ZSCORE.join(" ")} -prefix #{TTEST.sub("+tlrc.HEAD","")}")
end
CLOBBER.push(TTEST)
CLOBBER.push(TTEST.sub(".HEAD",".BRIK"))
CLOBBER.push(TTEST.sub(".HEAD",".BRIK.gz"))

file TTEST_ABS=> ZSCORE_ABS+[GROUP_MASK] do
  sh("3dttest++ -mask #{GROUP_MASK} -setA #{ZSCORE_ABS.join(" ")} -prefix #{TTEST_ABS.sub("+tlrc.HEAD","")}")
end
CLOBBER.push(TTEST_ABS)
CLOBBER.push(TTEST_ABS.sub(".HEAD",".BRIK"))
CLOBBER.push(TTEST_ABS.sub(".HEAD",".BRIK.gz"))

file TTEST_IGRAPH=> ZSCORE_IGRAPH+[GROUP_MASK] do
  sh("3dttest++ -mask #{GROUP_MASK} -setA #{ZSCORE_IGRAPH.join(" ")} -prefix #{TTEST_IGRAPH.sub("+tlrc.HEAD","")}")
end
CLOBBER.push(TTEST_IGRAPH)
CLOBBER.push(TTEST_IGRAPH.sub(".HEAD",".BRIK"))
CLOBBER.push(TTEST_IGRAPH.sub(".HEAD",".BRIK.gz"))

file TTESTBLUR => ZSCOREBLUR+[GROUP_MASK] do
  sh("3dttest++ -mask #{GROUP_MASK} -setA #{ZSCOREBLUR.join(" ")} -prefix #{TTESTBLUR.sub("+tlrc.HEAD","")}")
end
CLOBBER.push(TTESTBLUR)
CLOBBER.push(TTESTBLUR.sub(".HEAD",".BRIK"))
CLOBBER.push(TTESTBLUR.sub(".HEAD",".BRIK.gz"))

file TTESTBLUR_ABS => ZSCOREBLUR_ABS+[GROUP_MASK] do
  sh("3dttest++ -mask #{GROUP_MASK} -setA #{ZSCOREBLUR_ABS.join(" ")} -prefix #{TTESTBLUR_ABS.sub("+tlrc.HEAD","")}")
end
CLOBBER.push(TTESTBLUR_ABS)
CLOBBER.push(TTESTBLUR_ABS.sub(".HEAD",".BRIK"))
CLOBBER.push(TTESTBLUR_ABS.sub(".HEAD",".BRIK.gz"))

file TTESTBLUR_IGRAPH => ZSCOREBLUR_IGRAPH+[GROUP_MASK] do
  sh("3dttest++ -mask #{GROUP_MASK} -setA #{ZSCOREBLUR_IGRAPH.join(" ")} -prefix #{TTESTBLUR_IGRAPH.sub("+tlrc.HEAD","")}")
end
CLOBBER.push(TTESTBLUR_IGRAPH)
CLOBBER.push(TTESTBLUR_IGRAPH.sub(".HEAD",".BRIK"))
CLOBBER.push(TTESTBLUR_IGRAPH.sub(".HEAD",".BRIK.gz"))

file OVERLAP => AFNI do
  sh("3dmerge -1blur_fwhm 4.0 -gcount -prefix #{OVERLAP.sub("+tlrc.HEAD","")} #{AFNI.join(" ")}")
end
CLOBBER.push(OVERLAP)
CLOBBER.push(OVERLAP.sub(".HEAD",".BRIK"))
CLOBBER.push(OVERLAP.sub(".HEAD",".BRIK.gz"))

file PERM_MASK => [AVG,PERM] do
  tmpList = []
  (0..99).each do |iPerm|
    sh("3dcalc -a #{AVG} -b #{PERM}'[#{iPerm}]' -expr 'ispositive(a-b)' -prefix tmp#{iPerm}")
    tmpList.push("tmp#{iPerm}+tlrc.HEAD")
  end
  sh("3dmerge -gcount -prefix #{PERM_MASK.sub("+tlrc.HEAD","")} #{tmpList.join(" ")}")
  (0..99).each do |iPerm|
    rm("tmp#{iPerm}+tlrc.HEAD")
    rm("tmp#{iPerm}+tlrc.BRIK.gz") if File.exists?("tmp#{iPerm}+tlrc.BRIK.gz")
    rm("tmp#{iPerm}+tlrc.BRIK") if File.exists?("tmp#{iPerm}+tlrc.BRIK")
  end
end
CLOBBER.push(PERM_MASK)
CLOBBER.push(PERM_MASK.sub(".HEAD",".BRIK"))
CLOBBER.push(PERM_MASK.sub(".HEAD",".BRIK.gz"))

task :afni => AFNI
task :mean => AVG
task :overlap => OVERLAP
task :zscore => ZSCORE
task :ttest => [TTEST,TTESTBLUR]
task :ttest_abs => [TTEST_ABS,TTESTBLUR_ABS]
task :ttest_igraph => [TTEST_IGRAPH,TTESTBLUR_IGRAPH]
task :threshold => PERM_MASK
task :default => :threshold
task :rank => RANK + RANKBLUR
task :figures => PNG_TTEST + PNG_TTESTBLUR + PNG_NONPARAMETRIC + PNG_NONPARAMETRIC_MEAN
