#!/usr/bin/perl
package q_remote;
use strict;
use warnings;

############################################################
# declare variables
#-----------------------------------------------------------
use vars qw($scriptDir %config $cgi $params $isServer);
our ($nProjects, @projectNames, %project, @projectInputs);
our $projectsFile = "";
############################################################

############################################################
# get projects from file
#-----------------------------------------------------------
sub loadProjects {
    my $projectsDir = "$scriptDir/projects";
    unless(-d $projectsDir){
        mkdir $projectsDir;
        $isServer and qx|chmod -R ug+rw $projectsDir|;
    }
    $projectsFile = "$projectsDir/$config{user}.projects"; 
    my $projectsData = slurpFile($projectsFile);
    $projectsData =~ s/\r//g;
    extractProjects(\$projectsData, \my%projects);    
    setProjectNumber(\%projects); 
    foreach my $project(keys %projects){ 
        loadConfigLines(\$projects{$project}, \%{$config{projects}{$project}}) 
    }          
}
sub extractProjects {  # pull <project></project> blocks from q_remote.conf
    my ($projectsData, $projects) = @_;
    while($$projectsData =~ s|<project\s+(\w+?)>(.+?)</project>||s){ $$projects{$1} = $2 }
}
sub setProjectNumber {  # make sure at least one project is present
    my ($projects) = @_;
    $nProjects = scalar(keys %$projects);
    @projectNames = sort {$a cmp $b} keys %$projects;
    $nProjects > 1 and unshift @projectNames, "";
}
############################################################

############################################################
# project input and handling via web form
#-----------------------------------------------------------
sub setProject {  # on each page load, determine which project is in use
    %project = ();
    @projectInputs = ();
    createNewProject();
    deleteProject();
    $nProjects == 1 and $$params{project} = $projectNames[0];
    my $projectHash = $$params{project} ? $config{projects}{$$params{project}} : undef;
    %project = $projectHash ? %$projectHash : ();
    @projectInputs = projectInputs();  
    if($$params{projectChanged}){
        $$params{masterClass} = "";
        $$params{masterName} = "";
    }
} 
sub projectInputs {  # drop-down for selecting one of multiple projects
    my $padding = $$params{pageType} eq 'queue' ? 12 : 5;
    $padding = '&nbsp;' x $padding;
    my $new = "";
    $$params{pageType} eq 'queue' and
    $new = "&nbsp;&nbsp;<a href=javascript:newProject() title=\"Add a project directory for use in q remote\">Add</a>".
                         "\n<input type=hidden name=newProjectName id=newProjectName>\n".
	                     "\n<input type=hidden name=newProjectPath id=newProjectPath>\n".
	                     "&nbsp;&nbsp;<a href=javascript:deleteProject() title=\"Remove project from q remote (does NOT delete from server)\">Remove</a>".
	                     "\n<input type=hidden name=deleteProject id=deleteProject>\n";
    tablewrap2(["project$padding",
                 $cgi->popup_menu(-name=>"project", 
                                  -id=>"project",
                                  -values=>\@projectNames,
                                  -title=>"Select the project to use",
                                  -onChange=>"updateField('projectChanged')",
                                  -class=>"standard").
                 "\n<input type=hidden name=projectChanged id=projectChanged >\n" ,
                 $new] );                          
}
#-----------------------------------------------------------
sub createNewProject {  # create a new user project on server
    ($$params{newProjectName} and $$params{newProjectPath}) or return;
    getCurrentProjects(\my%projects);
    $projects{$$params{newProjectName}} = "\n    directory $$params{newProjectPath}\n";
    printProjects(\%projects);
    $$params{project} = $$params{newProjectName};  
}
sub deleteProject {  # delete a user project from the server
    $$params{deleteProject} or return;
    getCurrentProjects(\my%projects);
    delete $projects{$$params{project}};
    printProjects(\%projects);
    $$params{project} = "";
}
sub getCurrentProjects {
    my ($projects) = @_;
    my $projectsData = slurpFile($projectsFile);    
    extractProjects(\$projectsData, $projects);
}
sub printProjects {
    my ($projects) = @_;
    open my $outH, ">", $projectsFile or die "createNewProject error: could not open $projectsFile for writing\n";
    foreach my $project(keys %{$projects}){ 
        print $outH "<project $project>$$projects{$project}</project>\n";
    }
    close $outH;
    $isServer and qx|chmod ug+rw $projectsFile|;
    loadProjects(); 
}
############################################################

1;

