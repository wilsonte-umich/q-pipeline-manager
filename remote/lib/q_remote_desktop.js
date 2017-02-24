function setSizes(firstLoad){
    sshReturn = document.getElementById('sshReturn');    
    x = window.innerWidth - 455 - 10 ;
    y = window.innerHeight - 10 ;  
    sshReturn.style.width = x + 'px'; 
    sshReturn.style.height = y + 'px'; 
    sshWidth = x;
    sshHeight = y;     
    editBox = document.getElementById('editBox');
    if(editBox != null){
        offset = 10;   
        x -= offset;
        y -= offset;
        editBox.style.width = x + 'px'; 
        editBox.style.height = y + 'px';    
    }
}
