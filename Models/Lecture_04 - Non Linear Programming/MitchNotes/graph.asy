size(200,200);

pair[] verts = {
    (0,0), //A
    (-1,-1), //B
    (0,-1), //C
    (1,-1), //D
    (0,-2), //E
    (-1,-3), //F
    (0,-3), //G
    (1,-3), //H
    (1,-4), //I
    (1,-5), //J
    (0,-4), //K
    (1,-6), //L
    (1,-7), //M
    (1,-8), //N
    (2,-8), //O
    (3,-8), //P
    (4,-8) //Q
};

string[] names = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q"};

pair[] arcs = {
    (0,1),(0,2),(0,3),
    (1,4),
    (2,4),
    (3,4),
    (4,5),(4,6),(4,7),
    (5,13),
    (6,10),
    (6,8),
    (7,8),
    (8,9),
    (9,11),
    (10,12),
    (11,12),
    (12,13),
    (13,14),
    (14,15),
    (15,16)
};


for(int i = 0;i<names.length;++i){
    label(names[i],verts[i]);
}

for(pair p : arcs){
    pair source = verts[(int) p.x];
    pair dest = verts[(int) p.y];

    draw((1,1) -- (-5,-5));

    //draw(source--dest); 
}