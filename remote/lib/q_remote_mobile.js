function setSizes(firstLoad){
    if(firstLoad == 0){ return }  
    x = window.innerWidth;
    y = window.innerHeight;    
    editBox = document.getElementById('editBox');
    if(editBox != null){
        offset = 45;   
        x -= offset;
        y -= offset;
        editBox.style.width = x + 'px'; 
        editBox.style.height = y + 'px';    
    } 
}
