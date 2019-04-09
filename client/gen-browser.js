var app = new Vue({
    el: '#genom',


    data: {
        id: null,
        index: 0,
        alignment: null,
        alignment_count: 0,
        cigar: '',
        name: '',
        pos: 0,
        qualities: '',
        sequence: '',
        flags: '',
        atStart: true,
        atEnd: false
    },

    created: function () {
        this.count();
        this.read(1);
        this.id = 0;
    },

    methods: {
        count: function () {
            axios
                .get('http://localhost:8000/api/v1/count')
                .then(response => {
                    console.log(response.data);
                    this.alignment_count = Number(response.data);
                    console.log(this.alignment_count);
                });
        },
        read: function (index) {
            axios
                .get('http://localhost:8000/api/v1/alignment/' + index)
                .then(response => {
                console.log(response.data);

            this.name = response.data.name;
            this.pos = response.data.pos;
            this.index = response.data.index;
            this.cigar = response.data.cigar;
            this.sequence = response.data.sequence;
            this.qualities = response.data.qualities;

            var f = '';
            for (var flag in response.data.flags) {
                f += flag + ': ' + response.data.flags[flag] + ', ';
            }

            this.flags = f;

        })


        },

        first: function () {
            this.id = 0;
            this.read(1);
            this.atStart = true;
            this.atEnd = false;
        },

        next: function () {
            this.id += 1;
            this.read(this.id);
            this.atStart = false;
            this.atEnd = (this.id == this.alignment_count - 1)
        },

        prev: function () {
            this.id -= 1;
            this.read(this.id);
            this.atStart = (this.id == 0)
            this.atEnd = false;
        },

        last: function () {
            this.id = this.alignment_count - 1;
            this.read(this.id);
            this.atStart = false;
            this.atEnd = true;
        },
        clear: function () {
            this.name = '';
            this.pos = 0;
            this.index = 0;
            this.cigar = '';
            this.sequence = '';
            this.qualities = '';
            this.flags = '';
            this.id = 0;
            this.atStart = true;
        }
    }
})


