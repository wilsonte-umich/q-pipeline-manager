sshWidth = 100;
sshHeight = 100;
unselectedColor = '#DDDDDE';  // light grey
unselectedBorder = '';
selectedColor = '#EEEECC';  // lighter grey, touch of yellow
selectedBorder = 'inset';  
pageButtons = ['queueButton', 'shellButton'];
taskButtons = ['_editButton', '_submitButton', '_monitorButton', '_protectButton', '_publishButton',
               'headButton', 'tailButton', 'catButton', 'statButton', 'EditButton'];
subTaskButtons = ['masterButton', 'templateButton', 'helpButton',
                  'submitButton', 'extendButton', 'lockButton', 'rollbackButton', 'purgeButton',
                  'statusButton', 'clearButton', 'deleteButton', 'resubmitButton',
                  'protectButton', 'unprotectButton', 'backupButton', 'restoreButton',
                  'publishButton'];
isClicked = false;
clickedRow = false;
singleQuote = "__SINGLEQUOTE__";
newline = "__NEWLINE__";
function setPage(){
    setSizes(1);
    setPageButtons();
    setTaskButtons();
    setSubTaskButtons();
}
function setPageButtons() {
    pageType = document.getElementById('pageType').value;
    currentPageButtonID = pageType + "Button";
    currentPageButton = document.getElementById(currentPageButtonID);
    unselectButtons(pageButtons);
    selectButton(currentPageButton);
}
function setTaskButtons() {
    currentTask = document.getElementById('currentTask').value;
    currentTaskButtonID = currentTask + "Button";
    currentTaskButton = document.getElementById(currentTaskButtonID);
    unselectButtons(taskButtons);
    selectButton(currentTaskButton);
}
function setSubTaskButtons() {
    currentSubTask = document.getElementById('currentSubTask').value;
    currentSubTaskButtonID = currentSubTask + "Button";
    currentSubTaskButton = document.getElementById(currentSubTaskButtonID);
    unselectButtons(subTaskButtons);
    selectButton(currentSubTaskButton);
}
function unselectButtons(buttons) {
    for (i = 0; i < buttons.length; i++) {
        button = document.getElementById(buttons[i]);
        if(button != null){
            button.style.backgroundColor = unselectedColor;
            button.style.borderStyle = unselectedBorder;
        }
    }
}
function selectButton(button) {
    if(button != null){
        button.style.backgroundColor = selectedColor;
        button.style.borderStyle = selectedBorder;
    }
}
function resetPage(){
    document.getElementById('reset').value = 1; 
    document.forms['q_remote'].submit();
}
function setCacheDirection(direction){
    document.getElementById('cacheDirection').value = direction; 
    document.forms['q_remote'].submit();
}
function updateField(fieldID){
    document.getElementById(fieldID).value = 1; 
    document.forms['q_remote'].submit();  
}
function setPageType(pageType){
    currentPageType = document.getElementById('pageType');
    if(currentPageType.value == pageType){ return }
    currentPageType.value = pageType; 
    document.getElementById('pageTypeChanged').value = 1;
    setPageButtons();
    document.forms['q_remote'].submit();  
}
function setCurrentTask(task){
    currentTask = document.getElementById('currentTask');
    if(currentTask.value == task){ return }
    currentTask.value = task; 
    setTaskButtons();
    document.forms['q_remote'].submit();  
}
function setCurrentSubTask(subTask){
    currentSubTask = document.getElementById('currentSubTask');
    if(currentSubTask.value == subTask){ return }
    currentSubTask.value = subTask; 
    setSubTaskButtons();
    document.forms['q_remote'].submit();  
}
function executeTask(isDryRun,task,subTask){
    document.getElementById('executing').value = 1;
    if(isDryRun){ document.getElementById('isDryRun').value = 1 } 
    pageType = document.getElementById('pageType').value;
    if(pageType == 'queue'){
        if(task == '_edit'){ maskEditBox() }
        if(subTask == 'lock'){ document.getElementById('toggleLock').value = 1 }         
    } else if (pageType == 'shell') {
        if(task == 'Edit'){ maskEditBox() }
        shellCommand = document.getElementById('shellCommand');
        if(shellCommand != null){
            document.getElementById('maskedShellCommand').value = shellCommand.value.replace(/\n/g,' ; ');
        }
    }
    document.forms['q_remote'].submit();
}
function maskEditBox(){
    editBox = document.getElementById('editBox');
    if(editBox != null){
        str = editBox.value.replace(/'/g,singleQuote);
        document.getElementById('maskedEditBox').value = str.replace(/\n/g,newline);
    }
}
function getReport(jobID){
    document.getElementById('getReport').value = jobID;
    executeTask();
}
function changeRowColor(tableRow,highLight,rowID){
    if(rowID == 'deadRow') { return }
    isClickedID = "isClicked" + rowID;
    isClickedElement = document.getElementById(isClickedID);
    rowIsClicked = isClickedElement.value;
    if (highLight){
        tableRow.style.cursor='pointer';
        tableRow.style.backgroundColor = selectedColor; 
    } else {
        if(rowIsClicked == 1){
            tableRow.style.backgroundColor = unselectedColor;         
        } else {
            tableRow.style.backgroundColor = 'white';         
        }
    }
}
function processStatusListClick(event,tableCell,jobID){
    if(jobID == 'deadRow') { return }
    tableRow = tableCell.parentElement;
    isClickedID = "isClicked" + jobID;
    isClickedElement = document.getElementById(isClickedID);
    jobElement = document.getElementById('job');
    jobIDcomma = jobID + ",";    
    if(event.ctrlKey){
        tableRow.style.backgroundColor = 'white';
        isClickedElement.value = 0;
        if(jobElement != null && jobElement.value){
            if(jobElement.value.match(jobIDcomma)){ 
                jobElement.value = jobElement.value.replace(jobIDcomma,"") 
            }
        }
    } else {
        tableRow.style.backgroundColor = unselectedColor;
        isClickedElement.value = 1;
        if(jobElement != null){
            if(jobElement.value){
                if(!(jobElement.value.match(jobIDcomma))){ jobElement.value = jobElement.value + jobIDcomma; }
            } else {
                jobElement.value = jobIDcomma;
            }
        }
    }
}
function newMaster(masterType){
    newName = prompt("Name of the new master " + masterType);
    if (newName != null){
        flagID = "newMaster" + masterType;
        document.getElementById(flagID).value = newName;
        document.forms['q_remote'].submit(); 
    }
}
function newProject(){
    newPath = prompt("Enter the full path to a server directory with a 'masters' subdirectory");
    if (newPath != null){
        newName = prompt("Enter a project name for " + newPath, newPath);
        if (newName != null){
            document.getElementById('newProjectPath').value = newPath;        
            document.getElementById('newProjectName').value = newName;
            document.getElementById('projectChanged').value = 1;
            document.forms['q_remote'].submit(); 
        }
    }
}
function deleteProject(){
    response = confirm("Remove project entry'" + document.getElementById('project').value + "' ?");
    if (response == true){
        document.getElementById('deleteProject').value = 1;        
        document.forms['q_remote'].submit(); 
    }
}
function changeBrowseDir(){
    document.getElementById('browseDirChanged').value = 1;
}
