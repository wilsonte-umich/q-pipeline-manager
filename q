#!/usr/bin/perl
use strict;
use warnings;
$ENV{PATH} = "/home/wilsonte_lab/club_house/bin/q_versions/q-1.0.5/utilities/bedutil-1.0.0:/home/wilsonte_lab/club_house/bin/q_versions/q-1.0.5/utilities/slurp-1.0.0:/home/wilsonte_lab/club_house/bin/q_versions/q-1.0.5/utilities/stat-1.0.0:$ENV{PATH}";
our $version = '1.0.5';
our $perlPath = '/usr/bin/perl';
our $bashPath = '/bin/bash';
our $timePath = '/usr/bin/time';
our $timeVersion = '1.7';
our $memoryCorrection = '4';
our $memoryMessage = 'q: !! maxvmem value above is 4-fold too high due to known bug in GNU time utility !!';
our $qDir = '/home/wilsonte_lab/club_house/bin/q_versions/q-1.0.5';
our $libDir = '/home/wilsonte_lab/club_house/bin/q_versions/q-1.0.5/lib';
our $modulesDir = '/home/wilsonte_lab/club_house/bin/q_versions/q-1.0.5/modules';
our $qType = 'SGE';
our $schedulerDir = '/usr/local/gridengine/bin/linux-x64';
$ENV{Q_Q_TYPE} = 'SGE';
$ENV{Q_UTIL_DIR} = '/home/wilsonte_lab/club_house/bin/q_versions/q-1.0.5/utilities';
$ENV{Q_MOD_DIR} = '/home/wilsonte_lab/club_house/bin/q_versions/q-1.0.5/modules';
our %shellCommands = map { $_ => 1 } qw(pamtohtmltbl tcumttest pamcomp m4 gpgsplit dbilogstrip smartctl pitchplay ppmtoilbm pnmdepth gdk-pixbuf-query-loaders-64 doxygen ocs consoletype tclsh8.5 pack200 bigWigToBedGraph tzdata-update tload pammixinterlace TRFResult.pm start_udev automake lpr.cups vmstat lpinfo pnmpsnr h5repart editdiff chacl EmCommonCmdDriver.pm kar.2 eu-readelf runlevel newaliases.postfix pnmsmooth sra-dbcc.2 pnmconvol libvdb_jni.so sm-notify umount oraxml tcatest orabase ppdpo tinythread.cpp protoize h2ph sqlplus cracklib-format grep allneeded id qrsh ppmquant pnmscalefixed extproc swig neotoppm jmap idn hg19ToMm9.over.chain srvctl grpconv xfs_metadump fastjar djpeg libtool genhomedircon jbigtopnm gpg2 pmap papi_command_line thumbnail efibootmgr ppmbrighten addftinfo nenctool xfs_db glusterfsd unix2dos rsyslogd rgb2ycbcr pnmtojpeg msgunfmt grotty sldtoppm pgmtexture postsuper pi1toppm kpsexpand hg19ToHg18.over.chain rpc.nfsd plink rlatopam curl mkfs luac ppmtoppm pnmnlfilt otfview portrelease tiff2bw top db_recover zipinfo pbmtoxbm xzdec gtdownload gpasswd bowtie-inspect pdftotext MAKEDEV bcftools lesspipe.sh bonobo-activation-run-query zenity stapsh info readprofile tiff2rgba SeqDBI.pm lzgrep sgpio postcat vdb-dump.2 run_init pgmtopbm dbus-binding-tool grepdiff mergeBed uptime psidtopgm /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/blastdbcheck expandCols wnck-urgency-monitor mkfs.ext2 mbni-release-node pamfile vpddecode netreport pnmquant foomatic-fix-xml pamseq lsmod tchmttest rpmquery perldoc pnmtopng bedGraphToBigWig pambackground pacbio-load.2 ftp tnsping rpmdb timed-read ipcs grog pbmtomda extconv automount cpanp ncurses5-config ex pbmtoatk test-sra.2.3.2 gsl-config zipsplit dbfs_client gettext /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/convert2blastmask Rscript rename pambayer sgitopnm svnadmin deploync rmmod mksock vcfutils.pl eu-unstrip db_deadlock svnlook gs-906-linux_x86_64 sff-dump w gst-typefind install-info DFAMRecord.pm papi_clockres watch ppmwheel orion procob POST genorasdksh gif2tiff hexdump flipdiff gemtopnm fsck thinkjettopbm pear pacbio-load.2.3.2 sar od fsck.ext4 pnmhisteq urlgrabber ppmfade mdadm pamgradient certutil scrollkeeper-get-content-list dumphint slopBed papi_avail grub-md5-crypt fold dg4pwd xzcat kgmgr scrollkeeper-uninstall tfmtodit badblocks gsl-randist ldconfig xmlcatalog applygnupgdefaults pbmmake rarian-sk-preinstall qacct bc pgmdeshadow getopt isoinfo SecureUtil.pm makeindex blkdiscard extractKmers gffread cmuwmtopbm lzmainfo tgatoppm stap-merge kar autoheader pnmarith indxbib iconvconfig ck-launch-session sge_execd console-kit-daemon filterdiff pamchannel smime_keys intltool-extract foomatic-ppd-to-xml chown gserialver open_init_pty grmic dumpkeys dbus-cleanup-sockets rmid attr bamToFastq bowtie2-inspect infocmp EMSAConsoleCommon.pm abi-dump.2 rcsfreeze psfaddtable pnmtops pam_tally2 isodump pnmsplit gdb ppmtopuzz pamedge espdiff diagsetup nfsPatchPlugin.pm fgconsole foomatic-ppdfile pgmmedian grpunconv pnmnorm pbmtogo dircolors nenctool.2 dc aulast jarsigner ppmrough luseradd gctags x86_64-redhat-linux-gcc ppmtorgb3 kfod emdfail.command pdftex kfed fmt lid python2.6-config pstree.x11 mkfontscale pkexec showconsolefont sleep ccmake jrunscript pbmreduce papi_event_chooser mdmon blockdev ppmrainbow ps2pdf12 apropos nscd lpq.cups ppl-config dmraid getafm rsvg-convert p11-kit libpng12-config getsebool pal2rgb jstat rev ppmtopj papi_component_avail lsdiff grub setpci configuration-assistant.perl physe ls pinentry pbmtoepson fdformat gconf-merge-tree sprof pamtojpeg2k db_archive qrls fc-list gsdj500 msggrep cameratopam eqn2graph cpan LaunchEMagent.pm expdpO xfs_io lusermod mtvtoppm red gif2h5 grolj4 java-rmi.cgi map2gtf pbmtoybm mbni-mount objdump pbmto4425 wsgen pwd rcexplain.2.3.2 diff pam_timestamp_check podchecker crond groff kdbmeta h5unjam factor vdb-config.2 emutil.bat.template vimdiff whoami update-gdk-pixbuf-loaders sra-kar.2.3.2 findfs /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/segmasker x86_64-redhat-linux-c++ msgcomm eplain eu-elfcmp lsattr stap afs5log dvipdf newgrp perltex pam_console_apply kpsepath cjpeg df timeout xwdtopnm getOverlap pamthreshold foomatic-combo-xml iptc foomatic-configure route oldfind genisoimage foomatic-addpjloptions strace-log-merge msghack load_policy pkaction EMDeploy.pm texconfig-sys tcbmgr pamgauss dumpe2fs mkxauth lcsscan pear-0.9.10-bin-64 pampick Mail /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/legacy_blast.pl pnmmontage mktexpk env pbmtonokia gcov pgawk genhostid hpcdtoppm expand TRF.pm gcj-dbtool mkdir ifconfig msginit mbni-rsync rarian-sk-rebuild pfbtops emca xargs chsh iptables-multi-1.4.7 rarian-sk-config nfsiostat tiffcrop pamtopfm fsfreeze checkisomd5 podselect latex mktemp mii-diag scrollkeeper-extract addpart dbus-send brushtopbm mkfs.ext4dev gtbl EMAgent.pm pkill plipconfig grpck authconfig aleph ar tune2fs batch jar sra-sort.2.3.2 pbmupc otangle sff-dump.2 bed_to_juncs iptables-restore-1.4.7 segment_juncs pamsharpmap cifsiostat xzfgrep pamtosvg pftp jconsole peekfd unionBedGraphs json_pp wiggles unlzma pstopnm omega getent wrap tex rtsora echodo multiIntersectBed bowtie2-align library_stats ip6tables-save-1.4.7 unexpand open getfattr lxinst WUBlastSearchEngine.pm nameif nucBed dmevent_tool lpc.cups ssltap otfdump vdb-decrypt.2.3.2 ifenslave plymouthd runcon patgen illumina-load.2.3.2 rdisc gpic ppmrelief urlview ppmhist zcat svnsync tail gij fmtutil-sys osh h5mkgrp scrollkeeper-config db_upgrade ppmtopjxl qsub pwunconv epsffit jpegtran ck-log-system-restart unix_chkpwd emacs xjc tcttest glusterfs dbiprof gnomevfs-mkdir xbmtopbm /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/psiblast cvs rcs-checkin etags.emacs tracepath oerr git-shell cpp subtractBed sndfile-play floppy rhgb-client bwa DeCypherSearchEngine.pm who pnmcut unicode_start gdbus netmgr.ouibak rarian-sk-get-scripts grn javadoc ppdc pamfunc pbmtolps setmetamode pbmtoppa cg-load.2.3.2 gtnameserv aqxmlctl rcsdiff schemagen ip6tables-save gouldtoppm ps2ps pbmtoibm23xx hwloc-gather-topology cscope DateRepeats corelist jpegtopnm sgmlwhich lpmove bamToBed tie ppmtoterm udevd ArrayListIterator.pm qmake gdb-add-index vdb-copy.2.3.2 fastq-load.2.3.2 xzgrep dmp pbmclean vdb-validate.2.3.2 ssh bndlchk stty xml ausearch qrdel freetype-config poisson.TMP ppmcie readlink gnomevfs-cat jps closure_juncs logname /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/rpstblastn eu-addr2line relink pod2man ppmglobe iptables-save mail gorbd pnmtoplainpnm eu-ranlib xvminitoppm grmiregistry bedToBam cpio ppmtoicr arch date iptables-1.4.7 pamoil ppmtomitsu update-desktop-database busybox bzless setfiles delpart install lpc omfonts mkinitrd gmake stat groupBy sndfile-regtest ping xzcmp emacs-23.1 request-key tophat pbmtolj svn bzmore RepbaseRecord.pm dir ppdmerge h5debug iptables-xml-1.4.7 h5fc ca-legacy a2ping ul vimtutor appletviewer slattach combinediff pygtk-demo pdftosrc gtocheck gennttab imgtoppm rarian-example tic lockdev ppmnorm update-alternatives pbmtopsg3 lispmtopgm wsimport seq ppmtosixel bison exp sra_to_solid sln gs bzcmp kpsetool foomatic-rip mktexlsr switch_root xauth fastq-dump pgmramp domainname pktogf lxchknlb gettextize ps2ascii hwloc-assembler-remote mke2fs rdjpgcom /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/rpsblast pnmhistmap mkdumprd papi_mem_info grops modutil jvisualvm foomatic-kitload pjtoppm hpftodit firehose_get gnomevfs-ls scrollkeeper-gen-seriesid scriptreplay pnmcat h5dump pbmtomatrixorbital callgrind_annotate ip6tables tiff2pdf tangle sudo saslpasswd2 rmail pgmoil scrollkeeper-get-index-from-docpath unicode_stop HEAD zcmp ping6 sserver qrttoppm pkcheck zipgrep tiffsplit lchfn genoccish redhat_lsb_init ppdhtml gpg-agent psfstriptable instmodsh 411toppm xfs_estimate anytopnm newaliases wish8.5 ProcessRepeats sestatus mountstats pamsplit kpsewhere lzmadec pv.sh gnomevfs-monitor gftype statusnc prove kpseaccess cmake proxy kbd_mode racgeut mac2unix mag h5import kpsewhich sam_juncs pamarith bigWigToWig lprm.cups pamdepth rpmbuild gtk-update-icon-cache lstopo signtool qhost gnome-keyring qeseq fuser rcs cslatex kdump grmid ps2pdfwr gst-inspect-0.10 refseq-load.2 bzgrep rpm /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/tblastn pk2bm bibtex pgmnoise h5repack h5copy grubby ctangle pfbtopfa pdfopt pamtotiff extract_reads foomatic-preferred-driver qmod s2p doxytag clockdiff x86_64-redhat-linux-g++ plshprofO lbuilder lsof tcamgr SecureAgentCmds.pm /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/makeblastdb dos2unix isosize scrollkeeper-get-extended-content-list wftopfa cupstestppd saslauthd sendmail.postfix pbmtext pbmlife sha1sum h5cc rawtopgm unlink pbmpscale pbmtopgm manpath swapoff glookbib lacheck debugfs gconftool-2 ip6tables-1.4.7 gpgkey2ssh rpc.gssd split lprsetup.sh cpack ppmshadow base64 cacertdir_rehash losetup whatis pdfcrop xfs_info neqn nail dhclient iptunnel rpcgen pgmminkowski cracklib-packer illumina-load.2 reboot opl2ofm postlog snice bigBedToBed fipscheck abi-dump.2.3.2 lprm g++ hwloc-calc eu-stack pamditherbw rcexplain nosetests tcsh cupsfilter fstrim openssl texconfig-dialog qmon pod2usage ppmtolj zipnote fdisk pgrep ldd wget adapters pairToBed gensyslib commonenv rpc.mountd xpmtoppm atrm texconfig ppmdmkfont passmass ntpdate CompEMcore.pm dbhome tailf c89 ppmddumpfont cmsutil hwloc-distances cc cgquery fc-query fontinst lspci ssh-keygen impO rpc.svcgssd shred gst-xmllaunch bedToBigBed pamtotga comm tput git pnmpaste qtconfig-qt4 jstack pdftoppm tracepath6 mbni-sendthem pidof c2ph AgentLifeCycle.pm bonobo-activation-sysconf tiffmedian lastlog url_handler.sh ctest /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/blastdb_aliastool sff-load.2.3.2 vdb-dump.2.3.2 fsck.ext4dev stap-report ps4pdf byacc zipcloak ifdown expect dvips vim pdfcslatex /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/blastn atktopbm login grodvi auvirt flankBed sortBed sdiff numastat pbmtopk ps2pk mv eu-elflint lp.cups otflist git-upload-pack wc gst-feedback-0.10 gnomevfs-df otp2ocp ppmtopcx fastqc ppmtowinicon qconf iconvconfig.x86_64 bzcat sqlite3 mrftopbm chgrp SearchEngineI.pm fix_map_ordering accton iptables-multi setsid hwloc-distrib lzdiff ln ppmlabel NCBIBlastSearchEngine.pm DFAM.pm bam2fastx rediff bdftops pbmtowbmp adrciO pamstretch killall5 pwconv getenforce ppmtoyuvsplit dvi2fax otftobdf testrb gjarsigner vigr ppm2tiff tset pgmslice chrt catchsegv pbmtoptx gnome-text-editor portreserve wrjpgcom lzmore ott bzip2 as ctags sra-sort.2 gst-typefind-0.10 qdbusviewer du cg-load.2 etags modprobe vi pinentry-curses lmsgen sra-stat e2label ifnames pbmtog3 multiBamCov pamendian dracut gthumb-importer getkeycodes qualfa2fq.pl skill psbook mbni-chmod iconvconfig.i686 texsis xmlwf plurp colrm loadpsp ras2tiff mapUnique mf-nowin xzegrep gnome-open commonenv.template csscan lscpu RepeatMaskerConfig.tmpl spottopgm jsadebugd update-smart-drivedb ppmtompeg watchgnupg EMAgentPatch.pm javap mcookie pdftohtml rpm2cpio migspeed uuidgen pgmtopgm scp loadkeys bzip2recover lpr service tiffinfo pbmpage pamscale sysctl fiascotopnm checkmodule setup-nsssysinit.sh psresize scrollkeeper-get-toc-from-docpath ranlib ppmcolormask gennfgt netca selinuxenabled curl-config arping psnup foomatic-ppd-options basename gst-feedback a.out ddate postqueue phyzzx pamslice testsaslauthd icontopbm rarian-sk-extract init kdbmeta.2.3.2 winicontoppm finger lchsh gpgparsemail javafxpackager tbl modinfo wish update-ca-trust flock halt fstab-decode start matchpathcon lpasswd create-cracklib-dict setterm rcsclean pbmtopi3 blkid findmnt mdatopbm nenctool.2.3.2 pod2html logger msgen dirname sync sim_server pacbio-load autoupdate trcroute umount.nfs eu-make-debug-archive servertool numademo kgmgrO idlj helicos-load.2.3.2 unbuffer rsvg vmcore-dmesg unzipsfx bash foomatic-datafile mft papi_version gindxbib pamsharpness sra-dbcc semodule_expand patchAgtStPlugin.pm ppmtoarbtxt nohup pgmkernel eps2eps lpstat resize2fs mount.nfs metaflac msgfmt plshprof mount.nfs4 arp hwclock ntsysv weave maskFastaFromBed gcore mpost perlbug whereis genclntsh ppmtv pnmremap cupsdisable create_transcriptome_map scrollkeeper-preinstall ps2pdf dgmgrl sge_coshepherd reject setenforce foomatic-cleanupdrivers uname set_stat pamundice ghostscript nfsstat windowMaker sum rcs2log ipcalc papi_native_avail lslogins vdb-copy.2 tac lwp-download unprotoize ftp-rfc mutt cscope-indexer ip6tables-multi-1.4.7 mltex lambda scrollkeeper-get-cl /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/blastdbcmd nologin find ppmtoxpm rawtoppm printenv cvsbug cancel.cups mailq.postfix dbus-daemon setcap logsave pamstretch-gen luserdel man2html pbmtobbnbg python2 avcstat ncursesw5-config bedtools more rarian-sk-install vipw rapier tee lppasswd pnmtoddif cupsaccept h5diff clusterBed gtk-query-immodules-2.0-64 mkindex sh ccache-swig groups gemtopbm gsl-histogram semodule_link allcm scrollkeeper-rebuilddb pamlookup lpq StartAgent.pl authconfig-tui tcutest desktop-file-install look ngettext col qw pgmtolispm tabs SearchResultCollection.pm lsb_release bam_merge stop fast_mutex.h makedumpfile pcxtoppm mknod h2xs sra-dbcc.2.3.2 chcpu pamtopam tr dhclient-script mkfs.xfs makeinfo csplain fipshmac dmeventd flex++ musixflx piconv postconf libnetcfg infotopam yacc ctrlaltdel isodebug arpd hdifftopam pathchk captoinfo autoscan getkey wall postlock ucs2any gsbj bed12ToBed6 svnversion install-catalog oraxsl grep-changelog rpmsign macptopbm impdp canberra-gtk-play pamdice pbmtoescp2 mount.glusterfs rpc.statd xfs_freeze depmod rconfig db_printlog getfacl automake-1.11 dump-utmp pamtoxvmini ip6tables-restore pnmtotiffcmyk racgmain abi-load.2 dbus-uuidgen refer ifup koi8rxterm zip lookbib intltoolize ppmtojpeg tar kpartx foomatic-replaceoldprinterids status get_run_info mountpoint nice emutil sha256sum pfb2pfa ppmtoacad rpcinfo lgroupmod xfs_rtcp /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/blast_formatter updmap EmKeyCmds.pm emctl.template pltotf mesg gtf_to_sam rpmverify zegrep bowtie-build CompEMagent.pm R pre-grohtml packer intltool-update ppmtopict resizecons ipcrm ppmtomap vdb-encrypt.2 schema strip dbus-monitor pbmtodjvurle ps2ps2 qhold msgmerge mbni-chown RepeatProteinMask pamrgbatopng pcdovtoppm cg_diff ppmtoyuv get_stat dmidecode FastaDB.pm gm ppmtobmp dbshut installkernel [ lesskey cachefilesd a2p /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/blastx assistant_adp screen lzma rumakeindex eu-size segment locale selinuxconlist extcheck sudoreplay f2py sudoedit vdb-passwd.2.3.2 perlthanks clustalo-1.1.0-linux-64 foomatic-searchprinter msgconv lnstat readmult expr h5redeploy emacsclient tmpwatch autom4te vdb-passwd.2 closestBed gpg restorecon gst-xmlinspect enc2xs lgroupdel numactl lwp-rget sndfile-convert file ci less gsoelim genagtsh ppmpat recountdiff abi-dump scsi_id grub-crypt /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/blastp ebrowse git-receive-pack fc-match rtmon pnmtotiff pamtogif svndumpfilter loadjava last lwp-dump /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/windowmasker raw gendiff gst-inspect slogin loadunimap sysresv pamflip GET atq iptables-xml pgmmake bibtex8 xfs_mdrestore pamfixtrunc pk12util nm smartd umu usleep tiffset qresub Directory.pm x86_64 free pbmtomacp ssh-copy-id capsh smtp-sink compile_et slabtop make ddbugtopbm sndfile-metadata-get prefetch config_data rlog pbmtoepsi mailx gsnd uniq papi_xml_event_info rvi rman xfs_logprint ppmdim mkfs.ext4 cufflinks pgpring diskmon.bin pktype wrcO ppp-watch mount ident PubRef.pm gpgv2 pgmabel mbni-reserve-node pr iconv Path.pm contig_to_chr_coords gio-querymodules-64 AgentMisc.pm ip trcasst cat diffstat perlivp pnmtopclxl flex pgmtofs nfs_cache_getent mkfifo passwd pcregrep cbq configure pbmtoascii xfs_mkfile iptables-save-1.4.7 extjob swapon easy_install-2.6 db_dump extprocO updmap-sys foomatic-extract-text togglesebool umount.nfs4 pamaddnoise krb5-config kexec keyctl cancel RegisterTType.pm bzdiff gslp ppmtoleaf gff3ToGenePred reset timed-run refseq-load link /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/makembindex pf2afm jcontrol odvips pax vamp2 strings pic2graph pwdx plymouth-set-default-theme cupsreject md5sum font2c sra-pileup pdftops pamstack recode-sr-latin pcretest tcucodec plymouth expdp notify-send yuvtoppm javac getconf indent printafm lzegrep foomatic-getpjloptions isovfy rmic pamtohdiff mkofm tinythread.h ppdi cupsctl sAgentUtils.pm orbd mbni-du h5perf_serial update-mime-database cuffmerge aureport kill pamcut sirtopnm postfix pdf2ps etex sge_shadowd pdf2dsc rmiregistry RepbaseEMBL.pm colcrt qdbus vptovf fonttosfnt initctl visudo pdfetex sra-stat.2.3.2 dash gsdj pi3topbm lastb sputoppm hg18ToHg19.over.chain sqlldr fgrep pnmcrop rarian-sk-get-cl update-pciids cg_merge diff3 rcsmerge gst-launch gstack xxd netmgr mbni-receiveit zic mkisofs rarian-sk-migrate uuclient modulecmd tty ldattach setfont apt xzdiff vdb-config.2.3.2 glurp sndfile-info qdel /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/deltablast gencat callgrind_control ps2epsi epstopdf pampop9 qrstat pstops pamtodjvurle tnsping0 EMomsCmds.pm illumina-dump.2.3.2 amstex smtp-source sed pnmmargin dmstats echo gtoinfo ppmtoneo pl2pm /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/tblastx ac pnmtofiasco nsenter pango-view egrep AgentStatus.pm qd pairToPair lkbib ppmquantall faillock tnameserv e2fsck gnroff iptables pdfinfo mklost+found gneqn pnmcolormap svtools ed pydoc tiffdither pfmtopam ms_print pnmtorast python-config zeisstopnm wbmptopbm sushell setleds prefetch.2 ypdomainname kibitz tiffcmp ltrace weak-modules mbni-df pdffonts atd giftopnm papi_decode trcasst.ouibak hwloc-info fallocate soelim pamtopnm view sqlldrO ppmtopgm e2pall install_script ppmdist genomeCoverageBed nslcd tftopl kpsestat pluginviewer chcon allec touch lchage semodule pamstereogram expO cracklib-unpacker dbiproxy rvim owm zsoelim desktop-file-validate splain adrci uvProb bedToIgv nano c99 abi-load linux32 rlogin-cwd lpunlock pnmalias showchar lpoptions ps2pdf14 xsltproc addr2line toe i386 pslatex gnomevfs-rm partx proc native2ascii lpstat.cups gss-client write leaftoppm gftodvi xz groffer cweave wrc hwloc-ls securetty find2perl foomatic-printermap-to-gutenprint-xml scrollkeeper-install setfattr reload gprof fmtutil printf ppmcolors mbni-sequencing-core eu-nm valgrind-listener fsck.ext2 pgmnorm ionice escp2topbm pod2text gdbtui TRFSearchResult.pm gff2bed newer pnmtojbig pooltype find_single_isoforms ck-log-system-start scrollkeeper-update geqn nroff hwloc-assembler shuffleBed pgmedge gfortran rarian-sk-get-content-list localedef fusermount-glusterfs gzexe fixfiles extjobo telinit fixcvsdiff fstopgm slurp qping nfsidmap dmesg grub-terminfo pbmtomrf gthumb sndfile-cmp csslint-0.6 ld tagBam wigToBigWig pstree pnmtorle liftOver sndfile-metadata-set gkeytool sha512sum gjar cloog pnmindex imp rmdir kpsereadlink f2py.numpy tophat_reports sadf pnmscale fastqparser msgcmp bedutil nstat xfs_bmap EmctlCommon.pm at qalter hwloc-ps perl5.10.1 targetdeploy.pl gnsd bmp2tiff svnserve pwck cg-load cupsd pinky kbdrate cmp eyuvtoppm pamdeinterlace eu-findtextrel texhash geneziO scrollkeeper-get-toc-from-id DupMasker renice stap-prep bdf2gdfont.pl jdb splitdiff javaws vdb-copy pamenlarge gftopk semodule_deps ArrayList.pm tctmttest tc smooth mkhybrid csplit db_hotbackup pbmtozinc pbmmask orapki ppmmix true users igawk rarian-sk-update chpasswd ip6tables-multi grub-install implantisomd5 mm9ToHg19.over.chain gpg-connect-agent libtoolize outocp man rgb3toppm gst-xmlinspect-0.10 mkocp chvt annotateBed mllatex wipefs stdbuf chopt policytool pbmtocmuwm lxegen CompEMcentral.pm rview bgzip taskset mkhomedir_helper dd iostat pgmbentley gslj PRSearchResult.pm pbmto10x cpanp-run-perl postmap sfdisk tsort Matrix.pm foomatic-printjob pamtofits ppmspread prefetch.2.3.2 thumbpdf accept xml2-config yum pnmtosir ppmchange vdir bam-load.2.3.2 manweb mgrtopbm pc1toppm yes WUBlastXSearchEngine.pm qrsub semodule_package autoexpect jhat sg fc-cat restorecond nproc zdump memhog vdb-lock.2.3.2 rasttopnm xgettext gtf_juncs poisson pdfcsplain psfxtable ssh-add pnmstitch zfgrep ppmdither unix-lpr.sh sim_client uuserver pnmenlarge pgmmorphconv lgroupadd pnmcomp cuffcompare aulastlog ownership python dmsetup chfn aserver ppmflash ojvmtc ovp2ovf xzmore eu-strip windowBed tkprof pamtilt test msguniq lp assistant-qt4 shuf sa jv-convert pnminvert tcbtest column qq sendmail coraenv intltool-merge ip6tables-restore-1.4.7 mingetty insmod.static glib-compile-schemas rtacct atrun db_load gnomevfs-mv e2image infokey ethtool lpadmin ilbmtoppm sclient tzselect erb secon xfs_quota DESeq gst-launch-0.10 linksBed msgfilter cuffdiff hipstopgm samtools ofm2opl javah which jinfo runuser pstack grepjar pbmtogem c++ mkstore getcap fsck.cramfs checkpolicy xfs_admin serialver sbigtopgm dg4pwdO dnsdomainname ps db_dump185 biosdecode sasldblistusers2 pamtouil fc-cache lamed chage mkfontdir sam-dump.2.3.2 cumPoisson3 EMconnectorCmds.pm lzless perl sort texlinks sulogin EMDiag.pm postmulti utmpdump mktexmf makempx mkfs.cramfs b2m setkeycodes ppmntsc mkpasswd sha384sum pnmfile lnewusers pamsistoaglyph tcbmttest zdiff papi_multiplex_cost autopoint dmraid.static lzfgrep pnmgamma f95 unshare pnmtofits ps2frag hostid psfgettable sgepasswd pgpewrap cut pbmtox10bm rm pgmtoppm agetty xkibitz cracklib-check mount.tmpfs sge_shepherd awk pnmflip gcc enchant bioradtopgm pango-querymodules-64 lex pod2latex pamsummcol clear zforce c++filt raw2tiff ssh-agent mapscrn tclsh ximtoppm kdbmeta.2 udevadm join dvitomp mf unwrapdiff qtcsh jcmd crontab lsinitrd lzcat xfs_repair pammasksharpen pbmminkowski flac activation-client msgattrib eusm gss-server sra-pileup.2 align-info.2.3.2 long_spanning_reads setfacl setarch refseq-load.2.3.2 hwloc-bind tcfmgr getpcaps foomatic-nonumericalids poweroff libpng-config pktopbm cfdisk tchtest palmtopnm psed newusers lua troff sra-pileup.2.3.2 gnomevfs-info rpc.idmapd whiptail bunzip2 chattr /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/makeprofiledb h5c++ ss shutdown msgcat gedit g3topbm pnmshear gettext.sh sff-dump.2.3.2 ovf2ovp postdrop start-statd resize raid-check lastcomm sftp rnano lwp-request grolbp grefer tcfmttest xfs_fsr split_illumina xterm iptables-restore ybmtopbm script easy_install sedispol ruby zgrep nl papi_cost infotocap gpgconf genePredToGtf gnome-keyring-daemon psselect exportfs funzip tkprofO loadpspO fsck.xfs rsvg-view msgexec complementBed abi-load.2.3.2 h52gif pbmtoicon co gnomevfs-copy SearchResult.pm csh alternatives rmail.postfix /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/dustmasker SimpleBatcher.pm usernetctl CrossmatchSearchEngine.pm pbmtoplot osdbagrp patch pnmtoxwd uxterm ck-log-system-stop time aclocal-1.11 qselect mbchk kar.2.3.2 nisdomainname weather insmod asciitopgm ck-history sha224sum ausyscall java migratepages e2undo ppmmake vgdb fastq-dump.2.3.2 qseq2fastq dehtmldiff chkconfig mktexfmt dprofpp dropjava pgmhist cupsaddsmb fax2ps showmount lwp-mirror mpto sedismod dvired trcroute0 lessecho tiffcp gtroff tabix.py intersectBed cal gpgv gtf_to_fasta pkg-config mpstat mktextfm killall readelf dislocate Taxonomy.pm new-kernel-pkg xzless xfs_copy logrotate cp unxz bridge ppmtogif genezi emctl.pl ctie clock tctmgr rarian-sk-get-extended-content-list mkdict unzip vdb-encrypt.2.3.2 hostname intltool-prepare gnuplot rpcbind symfind fastaFromBed juncs_db gtar xsubpp dbstart gnuplot-wx coverageBed ipcmk aclocal pnmpad pamtompfont spctoppm size filefrag pngtopnm bdftopcf makewhatis h5ls helicos-load.2 rtcwake rpcdebug tifftopnm setsysfont lzcmp dump-acct vftovp berkeley_db_svc eu-ar postalias ppmtoeyuv pdfimages ssh-keyscan gzip mailq mask_sam ppmforge pamtooctaveimg zmore objcopy bmptoppm interdiff jpeg2ktopam gsftopk mii-tool linkshlib ptx rletopnm qstat ppmdcfont lstopo-no-graphics ppmdraw keytool h5stat db_codegen pic fsck.ext3 yuvsplittoppm cupstestdsc bowtie genclntst showkey ether-wake ncomp tchmgr rcexplain.2 tiffdump rmanO python2.6 tabix envsubst tiff2ps sge_qmaster prep_reads cg_annotate gunzip strace sra-stat.2 pnmnoraw pnmtopnm foomatic-perl-data netstat pnmtosgi rubibtex dm_dso_reg_tool fastq-dump.2 fax2tiff xfs_check truncate mkfs.ext3 pamperspective sra-kar.2 uidrvci orajaxb h5jam bmptopnm ppm3d ppmshift ppmtouil xfs_growfs gawk sshd SecureOMSCmds.pm false qlogin hunspell db_stat run-parts ifcfg ck-list-sessions foomatic-compiledb head cksum fc-scan selinuxdefcon pivot_root openvt RepeatMasker post-grohtml addgnupghome pnmtopalm randomBed /home/wilsonte_lab/club_house/bin/ncbi-blast-2.2.28+/bin/update_blastdb.pl deallocvt ipmaddr shasum AgentSubAgent.pm chmod xmllint xfs_ncheck crlutil pbmtomgr pstruct enchant-lsmod emwd.pl qsh su impdpO eu-strings cupsenable NCBIBlastXSearchEngine.pm oraenv rftp ps2pdf13 pidstat q eqn fitstopnm LineHash.pm db_checkpoint setsebool gpg-error tcftest znew cpuspeed namei postkick bowtie2 autoreconf papi_error_codes gst-xmllaunch-0.10 pgmcrater cpan2dist pbmtextps e2freefrag srf-load.2.3.2 bashbug-64 db_verify pamx qquota zless lsblk HMMERSearchEngine.pm paste sys-unconfig pamsumm eu-objdump pnmrotate pnmtile pbmtoln03 rsync tunelp chroot linux64 rarian-sk-gen-uuid pngtopam signver jstatd ppmtopi1 merge update-gtk-immodules pgmenhance restart valgrind git-upload-archive pdflatex pnminterp mkswap gpg-zip ojvmjava autoconf);
$| = 1;
require "$qDir/lib/main.pl";
qMain();
