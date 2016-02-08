
function uploader(callback){
     function handleFiles(e) {
        // console.log(e,this)
        e.stopPropagation();
        e.preventDefault();
        // d3.select('#files div').remove();

        var drop = false;
        var files = e.target.files; // for browse
        if (! files) {
            var files = e.dataTransfer.files; // for drag/drop
            drop = true;
            console.log('drop', drop);
        }

        var names = [];
        for (var i = 0; i < files.length; i++) {
            names.push(files[i]);
        };

        // just deal with first file
        var r = new FileReader();
        r.onload = function(e){
            callback(e.target.result, names);
        }
        r.readAsText(files[0]);
    }


    // select files by browsing or ....
    // ======================================


    var browse_inputs = document.getElementsByClassName('browse-input');
    for (var i = 0; i < browse_inputs.length; i++) {
        if (browse_inputs[i].addEventListener) browse_inputs[i].addEventListener('change', handleFiles, false);
    };

    var browse_buttons = document.getElementsByClassName('browse-button');
    for (var i = 0; i < browse_buttons.length; i++) {
        browse_buttons[i].addEventListener('click', function(e,b,c){
            // e.target == this
            this.parentNode.getElementsByClassName('browse-input')[0].click();
            // browse_buttons[i].previousSibling.previousSibling.click()
        }, false);
    };

    // ... drag and drop
    // ======================================

    function handleDragover(e) {
        e.stopPropagation();
        e.preventDefault();
        e.dataTransfer.dropEffect = 'copy';
    }

    var drop_divs = document.getElementsByClassName('drop');
    for (var i = 0; i < drop_divs.length; i++) {
        if(drop_divs[i].addEventListener) {
            drop_divs[i].addEventListener('dragenter', handleDragover, false);
            drop_divs[i].addEventListener('dragover', handleDragover, false);
            drop_divs[i].addEventListener('drop', handleFiles, false);
        }
    };
}

   