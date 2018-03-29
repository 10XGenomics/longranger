package test

import (
        "testing"
        . "loupe/formats"
        "encoding/json"
)


func TestCompress1(t * testing.T) {

        w, _ := NewCompressedWriter("testc.loupe", 4096);
        c1, _ := w.WriteChunk([]byte("Hello World! THis is a test!!!"));
        c2, _ := w.WriteChunk([]byte("Wow. Lets see if we can write data! To a file!"));
        l1 := []LoupeSection{c1,c2}
        
        j1, _:= json.Marshal(l1);

        w.WriteHeader(j1);
        w.Close();





}

