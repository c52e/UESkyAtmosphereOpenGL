#include "Textures.h"

#include "ImageLoader.h"
#include "GLReloadableProgram.h"

Textures::Textures() {
    {
        white_.Create(GL_TEXTURE_2D);
        glTextureStorage2D(white_.id(), 1, GL_RGBA8, 1, 1);
        uint8_t data[]{ 0xff,0xff,0xff,0xff };
        glTextureSubImage2D(white_.id(), 0, 0, 0, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, data);
    }
    {
        normal_.Create(GL_TEXTURE_2D);
        glTextureStorage2D(normal_.id(), 1, GL_RGBA8, 1, 1);
        uint8_t data[]{ 0x7f,0x7f,0xff,0xff };
        glTextureSubImage2D(normal_.id(), 0, 0, 0, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, data);
    }
    int x, y, n;
    {
        auto data = stbi_load_16_unique("../data/BlueNoise/64_64/HDR_L_0.png", &x, &y, &n, 0);

        blue_noise_.Create(GL_TEXTURE_2D);
        glTextureStorage2D(blue_noise_.id(), 1, GL_R16, x, y);
        glTextureSubImage2D(blue_noise_.id(), 0, 0, 0, x, y, GL_RED, GL_UNSIGNED_SHORT, data.get());
    }
    {
        auto data = stbi_load_unique("../data/NASA/lroc_color_poles_2k.jpg", &x, &y, &n, 0);

        moon_albedo_.Create(GL_TEXTURE_2D);
        glTextureStorage2D(moon_albedo_.id(), GetMipmapLevels(x, y), GL_SRGB8, x, y);
        glTextureSubImage2D(moon_albedo_.id(), 0, 0, 0, x, y, GL_RGB, GL_UNSIGNED_BYTE, data.get());
    }
    glGenerateTextureMipmap(moon_albedo_.id());
    {
        auto data = stbi_load_unique("../data/NASA/moon_normal_2k.jpg", &x, &y, &n, 0);

        // The moon normal is stored with DirectX normal direction
        // Change to glTF (OpenGL) normal direction
        for (int i = 0; i < x * y; ++i)
            data.get()[i * 3 + 1] = 255 - data.get()[i * 3 + 1];

        moon_normal_.Create(GL_TEXTURE_2D);
        glTextureStorage2D(moon_normal_.id(), GetMipmapLevels(x, y), GL_RGB8, x, y);
        glTextureSubImage2D(moon_normal_.id(), 0, 0, 0, x, y, GL_RGB, GL_UNSIGNED_BYTE, data.get());
    }
    glGenerateTextureMipmap(moon_normal_.id());
    {
        auto data = stbi_load_unique("../data/NASA/starmap_2020_4k.jpg", &x, &y, &n, 0);

        star_luminance_.Create(GL_TEXTURE_2D);
        glTextureStorage2D(star_luminance_.id(), GetMipmapLevels(x, y), GL_SRGB8, x, y);
        glTextureSubImage2D(star_luminance_.id(), 0, 0, 0, x, y, GL_RGB, GL_UNSIGNED_BYTE, data.get());
    }
    glGenerateTextureMipmap(star_luminance_.id());
    {
        auto data = stbi_load_unique("../data/NASA/world.topo.bathy.200401.3x5400x2700.jpg", &x, &y, &n, 0);

        earth_albedo_.Create(GL_TEXTURE_2D);
        glTextureStorage2D(earth_albedo_.id(), GetMipmapLevels(x, y), GL_SRGB8, x, y);
        glTextureSubImage2D(earth_albedo_.id(), 0, 0, 0, x, y, GL_RGB, GL_UNSIGNED_BYTE, data.get());
    }
    glGenerateTextureMipmap(earth_albedo_.id());

    {
        constexpr int kSizeX = 512;
        constexpr int kSizeY = 512;
        env_brdf_lut_.Create(GL_TEXTURE_2D);
        glTextureStorage2D(env_brdf_lut_.id(), 1, GL_RG16, kSizeX, kSizeY);
        GLReloadableComputeProgram program = {
            "../shaders/Base/EnvBRDFLut.comp",
            {{8, 8}}
        };
        glUseProgram(program.id());
        GLBindImageTextures({ env_brdf_lut_.id() });
        program.Dispatch({ kSizeX, kSizeY });
        glMemoryBarrier(GL_TEXTURE_FETCH_BARRIER_BIT);
    }
}

