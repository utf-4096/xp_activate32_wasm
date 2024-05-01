export default class XPA {
    static Errors = {
        '1': 'Too short',
        '2': 'Too large',
        '3': 'Invalid character',
        '4': 'Invalid check digit',
        '5': 'Unknown version',
        '6': 'Unlucky',
    };

    /** @type {string} */
    #binary_path = null;

    /** @type {?WebAssembly.Instance} */
    #instance = null;

    /** @type {Number} */
    #pBuffer = -1;

    constructor(
        /* path to the xpa wasm module */
        binary = './assets/wasm/xpa.wasm',
    ) {
        this.#binary_path = binary;
    }

    #write_ntstring(index, text) {
        const buffer = new Uint8Array(this.#instance.exports.memory.buffer);
        const encoder = new TextEncoder();
        const data = encoder.encode(text);
        buffer.set(data, index);
        buffer[index + data.byteLength] = 0;

        return index + data.byteLength + 1;
    }

    #read_ntstring(index) {
        const buffer = new Uint8Array(this.#instance.exports.memory.buffer);
        const decoder = new TextDecoder();

        let read = [];
        for (let i = 0; buffer[index + i] !== 0; i++) {
            read[i] = buffer[index + i];
        }

        return decoder.decode(new Uint8Array(read));
    }

    async boot() {
        if (this.#instance) {
            throw new Error('xpa already booted');
        }

        let instance;
        if ('instantiateStreaming' in WebAssembly) {
            const stream = fetch(this.#binary_path);
            instance = (await WebAssembly.instantiateStreaming(stream)).instance;
        } else {
            const code = await fetch(this.#binary_path).then(r => r.arrayBuffer());
            instance = (await WebAssembly.instantiate(code)).instance;
        }

        this.#pBuffer = instance.exports.get_buffer();
        this.#instance = instance;
    }

    generate(key) {
        if (!this.#instance) {
            throw new Error('boot() needs to be called first');
        }

        /* we don't want to overflow the buffer, so throw immediately */
        if (key.length > 62) {
            throw new Error(XPA.Errors['2']);
        }

        const pInput = this.#pBuffer;
        const pInputEnd = this.#write_ntstring(pInput, key);
        const pOutput = pInputEnd + 1;

        const status = this.#instance.exports.generate(pInput, pOutput);
        if (status !== 0) {
            throw new Error(XPA.Errors[status] ?? `Unknown error (${status})`);
        }

        return this.#read_ntstring(pOutput);
    }
}
